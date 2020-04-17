/*
	Copyright 2020 Marco Carolla
	
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    */

#include "obliviousTransfer.c"



int main(int argc, int *argv[])
{

    int Alice_Secrets = 2;

    if(argc > 1)
    {
        Alice_Secrets = atoi(argv[1]);
    }

    printf("1-out-of-%d SIDH-based OT\n", Alice_Secrets);

    // Test Benchmark
    int total_BenchmarkCycles = 10;
    unsigned long long total_Cycles = 0, start_Cycles, end_Cycles;

    for(int i = 0; i < total_BenchmarkCycles; i++)
    {
        start_Cycles = cpucycles();


        f2elm_t A24plus = {0}, C24 = {0}, A_f2element = {0};

        // Initialize constants: A24plus = A_f2element+2C, C24 = 4C, where A_f2element=6, C=1
        fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
        fp2add(A24plus, A24plus, A24plus);
        fp2add(A24plus, A24plus, C24);
        fp2add(A24plus, C24, A_f2element);
        fp2add(C24, C24, A24plus);

        /* - - - Part 1 - - - Alice - - - */
        /* - - - Part 1 - - - Alice - - - */
        /* - - - Part 1 - - - Alice - - - */

        unsigned char* secrets[Alice_Secrets][MSG_BYTES];
        point_proj_t Alice_2n_Basis[3*Alice_Secrets], Alice_isogenyPoints[3*Alice_Secrets], Alice_Base_phiT[Alice_Secrets];
        f2elm_t Alice_isogenousCurves[2*Alice_Secrets]; //This is formatted as pairs (A24plus, C24) for each computed curve

        for(int i = 0; i < Alice_Secrets; i++)
        {
            //Trivial initialization of Alice's secrets
            //TODO: You may want to change this
            *secrets[i] = rand();
            
            generate_2n_entangledBasis(A_f2element, Alice_2n_Basis[3*i], Alice_2n_Basis[3*i +1]);

            pointDifference(
                    Alice_2n_Basis[3*i +1], Alice_2n_Basis[3*i],
                    A_f2element, A24plus, C24,
                    Alice_2n_Basis[3*i +2]
                    );
            //At this point we have:
            //R0 = Alice_2n_Basis[3*i], first 2^n base point
            //T0 = Alice_2n_Basis[3*i + 1], second 2^n base point
            //Alice_2n_Basis[3*i + 2], their difference T0 - R0

            get_2n_isogenousCurve_4pt_CurveCoeff(
                    Alice_2n_Basis[3*i], Alice_2n_Basis[3*i +1],
                    Alice_isogenyPoints[3*i], Alice_isogenyPoints[3*i +1], Alice_isogenyPoints[3*i +2], Alice_Base_phiT[i],
                    A24plus, C24,
                    Alice_isogenousCurves[2*i], Alice_isogenousCurves[2*i +1]
                    );
            //Now we have:
            //3 Alice_isogenyPoints[] indexed 3i, 3i+1, 3i+2
            //1 Alice_Base_phiT[] which will be used in ALICE PART 3
            //1 Alice_isogenousCurves[] with its coefficients at 2i and 2i+1
        }

        //Now Alice may either:
        //send Bob Alice_isogenousCurves[] + Alice_isogenyPoints[]
        //or compress the Alice_isogenyPoints[] and send them to Bob
        //
        //TODO: compressione punti come riportato in sidh.c, riga 295-300.
        // La compressione va fatta Alice_Secrets volte in modo da comprimere tutti punti calcolati in Alice_isogenyPoints[]


        /* - - - Part 2 - - - Bob - - - */
        /* - - - Part 2 - - - Bob - - - */
        /* - - - Part 2 - - - Bob - - - */

        int bob_ChosenSecret_k = rand() % Alice_Secrets; //This should be chosen by Bob, here it's random
        unsigned char* bobSecret_3m = malloc(SECRETKEY_B_BYTES);
        random_mod_order_B(bobSecret_3m); //Bob secret key

        f2elm_t xPB, xQB, xRB;
        point_proj_t isogenyKernel;

        init_basis((digit_t*)B_gen, xPB, xQB, xRB);
        LADDER3PT(xPB, xQB, xRB, (digit_t*)bobSecret_3m, BOB, isogenyKernel, A_f2element); //isogenyKernel = P + (bobSecret_3m * Q)

        
        f2elm_t isogen3m_A24plus, isogen3m_C24; //3m-isogenous curve coefficients
        get_3m_isogenousCurve_Only(isogenyKernel, A24plus, C24, isogen3m_A24plus, isogen3m_C24);

        f2elm_t Bob_j_Inv;
        j_inv(isogen3m_A24plus, isogen3m_C24, Bob_j_Inv); //Bob now has his j-inv

        point_proj_t Bob_2n_Basis[3*Alice_Secrets], Bob_isogenyPoints[3];
        point_full_proj_t Bob_2n_Full_basis[2];

        int weilPairing_Outputs = 2;
        f2elm_t weilPairing_Value[4] = {0}, Bob_f2elem_k = {0};
        f2elm_t weierstrass_a_coefficient, weierstrass_b_coefficient;

        fp2add(Alice_isogenousCurves[0], Alice_isogenousCurves[1], Bob_f2elem_k);

        Monty2Weier(Bob_f2elem_k, weierstrass_a_coefficient, weierstrass_b_coefficient);

        generate_2n_FULL_entangledBasis(Bob_f2elem_k, Bob_2n_Full_basis[0], Bob_2n_Full_basis[1]);
        Convert_Full_to_Short(Bob_2n_Full_basis[0], Bob_2n_Basis[0]);
        Convert_Full_to_Short(Bob_2n_Full_basis[1], Bob_2n_Basis[1]);

        compute_Tate_Order2n(Bob_2n_Full_basis[0], Bob_2n_Full_basis[1], 
                            weierstrass_a_coefficient, Bob_f2elem_k, 
                            weilPairing_Outputs, weilPairing_Value);


        //Generate Alice_Secrets bases all of order 2^n and all having the SAME WEIL PAIRING
        for(int i = 1; i < Alice_Secrets; i++)
        {
            fp2add(Alice_isogenousCurves[2*i], Alice_isogenousCurves[2*i +1], Bob_f2elem_k);

            do
            {    
                //TODO: This function should take "Bob_f2elem_k" instead of "A_f2element" but it runs in an infinite loop
                generate_2n_FULL_entangledBasis(A_f2element, Bob_2n_Full_basis[0], Bob_2n_Full_basis[1]);
                compute_Tate_Order2n(Bob_2n_Full_basis[0], Bob_2n_Full_basis[1], 
                                    weierstrass_a_coefficient, Bob_f2elem_k, 
                                    weilPairing_Outputs, weilPairing_Value+2);
            } 
            while (
                weilPairing_Value[2] == weilPairing_Value[0]
                &&
                weilPairing_Value[3] == weilPairing_Value[1]
                );
        
            Convert_Full_to_Short(Bob_2n_Full_basis[0], Bob_2n_Basis[3 * i]);
            Convert_Full_to_Short(Bob_2n_Full_basis[1], Bob_2n_Basis[3 * i + 1]);

            pointDifference(
                    Bob_2n_Basis[3 * i + 1], Bob_2n_Basis[3 * i],
                    Bob_f2elem_k, Alice_isogenousCurves[2*i], Alice_isogenousCurves[2*i +1],
                    Bob_2n_Basis[3 * i + 2]
                );
        }

        //Now we have to compute one kernel and isogenous curve for each Alice_Secrets
        
        //isogenyKernel = P + (bobSecret_3m * Q)
        LADDER3PT(Alice_isogenyPoints[3*bob_ChosenSecret_k], Alice_isogenyPoints[3*bob_ChosenSecret_k +1], Alice_isogenyPoints[3*bob_ChosenSecret_k +2], (digit_t*)bobSecret_3m, BOB, isogenyKernel, A_f2element); 

        get_2n_isogenousCurve_3pt(
                isogenyKernel,
                Bob_2n_Basis[3*bob_ChosenSecret_k], Bob_2n_Basis[3*bob_ChosenSecret_k +1], Bob_2n_Basis[3*bob_ChosenSecret_k +2], //Points to "isogenize"
                Bob_isogenyPoints[0], Bob_isogenyPoints[1], Bob_isogenyPoints[2], //Their isogeny
                A24plus, C24
        );
        //Now we have:
        //3 Bob_isogenyPoints[]: U, V and V-U
        //1 Bob_isogenousCurves[] with its coefficients at 0 and 1 but this is not used further on

        // Now Bob has to send Alice: Bob_2n_Basis[], Bob_isogenyPoints[]
        // Similarly to Alice's, Bob may instead send Bob_2n_Basis[] and a compressed version of Bob_isogenyPoints[]


        // /* - - - Part 3 - - - Alice - - - */
        // /* - - - Part 3 - - - Alice - - - */
        // /* - - - Part 3 - - - Alice - - - */


        //x, y linear combination such that Alice_Base_phiT = x*U + y*V
        // NOTE: This does NOT work. It just gives 2 random values to x and y. Check linearCombination_Random() to better understant what's going on
        unsigned char* x_y[2 * Alice_Secrets];
        linearCombination_Random(Alice_Secrets, Bob_2n_Basis, Alice_Base_phiT, A_f2element, A24plus, C24, x_y);

        point_proj_t xU, difference, final_Isogeny_Kernel;
        f2elm_t Alice_j_Inv[Alice_Secrets];
        unsigned char* encrypted_secrets[Alice_Secrets][MSG_BYTES];

        for(int i = 0; i < Alice_Secrets; i++)
        {
            pointMultiplication(x_y[2*i], Bob_isogenyPoints[0], A_f2element, A24plus, C24, xU);      
            pointDifference(Bob_isogenyPoints[1], xU, A_f2element, A24plus, C24, difference); //difference = V - x*U
            LADDER3PT(xU->X, Bob_isogenyPoints[1]->X, difference->X, 
                        (digit_t*)&x_y[2*i +1], ALICE, final_Isogeny_Kernel, A_f2element);  //final_Isogeny_Kernel = x*U + y*V

            //Here we start replacing OLD Alice_isogenousCurves with the new ones so we don't initialize a new array
            get_2n_isogenousCurve_Only(
                                    final_Isogeny_Kernel, 
                                    A24plus, C24, 
                                    Alice_isogenousCurves[2*i], Alice_isogenousCurves[2*i +1] //Get the isogenous curve
                                    ); 

            j_inv(Alice_isogenousCurves[2*i], Alice_isogenousCurves[2*i +1], Alice_j_Inv[i]); //Alice now has her j-inv

            Enc(secrets[i], Alice_j_Inv[i], encrypted_secrets[i]);
        }
        //Alice sends encrypted_secrets[] to Bob


        /* - - - Part 4 - - - Bob - - - */
        /* - - - Part 4 - - - Bob - - - */
        /* - - - Part 4 - - - Bob - - - */


        unsigned char* decrypted_secret[MSG_BYTES];
        Enc(encrypted_secrets[bob_ChosenSecret_k], Bob_j_Inv, decrypted_secret);

        //Quit benckmark test
        end_Cycles = cpucycles();
        total_Cycles = total_Cycles+(end_Cycles-start_Cycles);

        // printf("Cycle %d:   %10lld  cycles\n", i+1, end_Cycles-start_Cycles);
        printf("%10lld\n", end_Cycles-start_Cycles);
        
        int result = check_Result(secrets[bob_ChosenSecret_k], decrypted_secret, MSG_BYTES);
        printf("Bob ha calcolato lo stesso segreto di Alice? -> %d\n", result);
    }

    printf("\nProtocol runs in %10lld cycles\n", total_Cycles/total_BenchmarkCycles);
    return 0;
}