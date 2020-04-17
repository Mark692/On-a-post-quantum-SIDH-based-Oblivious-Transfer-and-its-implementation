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

    // Test Benchmark
    int total_BenchmarkCycles = 20;
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


        /* - - - Part 0 - - - Public Points - - - */
        /* - - - Part 0 - - - Public Points - - - */
        /* - - - Part 0 - - - Public Points - - - */
        
        //P_A, Q_A, and their difference
        point_proj_t P_A, Q_A, QA_minus_PA; // = Q_A - P_A
        generate_2n_entangledBasis(A_f2element, P_A, Q_A);
        pointDifference(P_A, Q_A,
                        A_f2element, A24plus, C24,
                        QA_minus_PA);

        //P_B, Q_B, and their difference
        point_full_proj_t P_B_full, Q_B_full;
        point_proj_t P_B, Q_B, QB_minus_PB;
        unsigned int rs[2] = {0};

        BuildOrdinaryE3nBasis(A_f2element, P_B_full, Q_B_full, rs);
        Convert_Full_to_Short(P_B_full, P_B);
        Convert_Full_to_Short(Q_B_full, Q_B);
        pointDifference(P_B, Q_B,
                        A_f2element, A24plus, C24,
                        QB_minus_PB);



        /* - - - Part 1 - - - Alice - - - */
        /* - - - Part 1 - - - Alice - - - */
        /* - - - Part 1 - - - Alice - - - */

        unsigned char* Alice_r_A = malloc(SECRETKEY_A_BYTES);
        random_mod_order_A(Alice_r_A); //Alice integer r_A

        point_proj_t phiP_A, phiQ_A, phi_QmP, //Images of P_A, Q_A, QA_minus_PA
                     PA_plus_rQA; // = P_A + Alice_r_A*Q_A;

        LADDER3PT(P_A, Q_A, QA_minus_PA, (digit_t*)Alice_r_A, ALICE, PA_plus_rQA, A_f2element);

        f2elm_t curveE_A_coefficients[2] = {0};
        get_2n_isogenousCurve_3pt_CurveCoeff(PA_plus_rQA, P_A, Q_A, QA_minus_PA, phiP_A, phiQ_A, phi_QmP, A24plus, C24, *curveE_A_coefficients[0], *curveE_A_coefficients[1]);

        //Alice public key is the triple {curveE_A_coefficients, phiP_A, phiQ_A}
        //TODO: compressione punti come riportato in sidh.c, riga 295-300.


        /* - - - Part 2 - - - Bob - - - */
        /* - - - Part 2 - - - Bob - - - */
        /* - - - Part 2 - - - Bob - - - */

        int sigma = rand() % Alice_Secrets; //This should be chosen by Bob, here it's random
        unsigned char* Bob_r_B = malloc(SECRETKEY_B_BYTES);
        random_mod_order_B(Bob_r_B); 

        point_proj_t phiP_B, phiQ_B, phiB_QmP, //Images of P_A, Q_A, QA_minus_PA
                     PB_plus_rQB; // = P_B + Bob_r_B*Q_B;

        LADDER3PT(P_B, Q_B, QB_minus_PB, (digit_t*)Bob_r_B, BOB, PB_plus_rQB, A_f2element);

        f2elm_t curveE_B_coefficients[2] = {0};
        get_3m_isogenousCurve_3pt_Coeff(PB_plus_rQB, P_B, Q_B, QB_minus_PB, A24plus, C24, phiP_B, phiQ_B, phiB_QmP, curveE_B_coefficients[0], curveE_B_coefficients[1]);

        //Bob's FIRST public key is the triple {curveE_B_coefficients, phiP_B, phiQ_B}

        point_full_proj_t U_full, V_full;
        point_proj_t U, V;
        f2elm_t curveE_B_coeff_A;
        unsigned int rs2[2] = {0};

        //TODO: This code should be debugged and then replace the BuildOrdinaryE3nBasis(A_f2element, U_full, V_full, rs2);
        //fp2add(curveE_B_coefficients[0], curveE_B_coefficients[1], curveE_B_coeff_A);
        // BuildOrdinaryE3nBasis(curveE_B_coeff_A, U_full, V_full, rs2);
        BuildOrdinaryE3nBasis(A_f2element, U_full, V_full, rs2);
        fp2copy(A_f2element, curveE_B_coeff_A);

        Convert_Full_to_Short(U_full, U);
        Convert_Full_to_Short(V_full, V);

        point_proj_t sigma_U, sigma_V, GB_tilde, HB_tilde;
        //Attento che sta funzione prende char e non int
        pointMultiplication(sigma, U, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], sigma_U);      
        pointMultiplication(sigma, V, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], sigma_V);      

        pointDifference(phiP_B, sigma_U, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], GB_tilde);
        pointDifference(phiQ_B, sigma_V, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], HB_tilde);

        //Bob's SECOND public key is the triple {curveE_B_coefficients, GB_tilde, HB_tilde}


        // /* - - - Part 3 - - - Alice - - - */
        // /* - - - Part 3 - - - Alice - - - */
        // /* - - - Part 3 - - - Alice - - - */

        unsigned char* secrets[Alice_Secrets][MSG_BYTES];
        unsigned char* encrypted_secrets[Alice_Secrets][MSG_BYTES];
        point_proj_t iU, iV, sum1, sum2, sumDiff, isogenyKernel;
        f2elm_t image_A24, image_C24, Alice_jinv;

        for(int i = 0; i < Alice_Secrets; i++)
        {
            //Trivial initialization of Alice's secrets
            //TODO: You may want to change this
            *secrets[i] = rand();

            pointMultiplication(i, U, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], iU);      
            pointMultiplication(i, V, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], iV);      

            pointAddition(GB_tilde, iU, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], sum1);
            pointAddition(HB_tilde, iV, curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1], sum2);

            pointDifference(sum1, sum2,
                            curveE_B_coeff_A, curveE_B_coefficients[0], curveE_B_coefficients[1],
                            sumDiff);

            LADDER3PT(sum1, sum2, sumDiff, (digit_t*)Alice_r_A, ALICE, isogenyKernel, curveE_B_coeff_A);
            get_3m_isogenousCurve_Only(isogenyKernel, curveE_B_coefficients[0], curveE_B_coefficients[1], image_A24, image_C24);

            j_inv(image_A24, image_C24, Alice_jinv);
            Enc(secrets[i], Alice_jinv, encrypted_secrets[i]);
        }
        //Alice sends encrypted_secrets[] to Bob


        /* - - - Part 4 - - - Bob - - - */
        /* - - - Part 4 - - - Bob - - - */
        /* - - - Part 4 - - - Bob - - - */

        unsigned char* decrypted_secret[MSG_BYTES];

        point_proj_t finalIsogeny;
        f2elm_t curveE_A_coeff_A, finalA24, final24, Bob_j_Inv;
        fp2add(curveE_A_coefficients[0], curveE_A_coefficients[1], curveE_A_coeff_A);
        LADDER3PT(phiP_A, phiQ_A, phi_QmP, (digit_t*)Bob_r_B, BOB, finalIsogeny, curveE_A_coeff_A);
        get_2n_isogenousCurve_Only(finalIsogeny, curveE_A_coefficients[0], curveE_A_coefficients[1], finalA24, final24);

        j_inv(finalA24, final24, Bob_j_Inv);
        Enc(encrypted_secrets[sigma], Bob_j_Inv, decrypted_secret);

        //Quit benckmark test
        end_Cycles = cpucycles();
        total_Cycles = total_Cycles+(end_Cycles-start_Cycles);

        // printf("Cycle %d:   %10lld  cycles\n", i+1, end_Cycles-start_Cycles);
        printf("Cicli CPU necessari: %10lld", end_Cycles-start_Cycles);
        
        
        int result = check_Result(secrets[sigma], decrypted_secret, MSG_BYTES);
        printf(", Bob ha decifrato correttamente il segreto? -> %d\n\n", result);
    }

    printf("\nMedia cicli necessari: %10lld cycles\n", total_Cycles/total_BenchmarkCycles);
    return 0;
}