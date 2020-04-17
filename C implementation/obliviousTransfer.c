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

   
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> //To printf("%" PRIx64 "\n", variabile_uint64_t);
#include "lib/P503/api.h"
#include "lib/config.h"
#include "lib/P503/P503_internal.h"
#include "lib/P503/P503.c"

#include "obliviousTransfer.h"


digit_t projective_Z_AtInfinity503[NWORDS_FIELD] = 
{
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xabffffffffffffff,
    0x13085bda2211e7a0, 0x1b9bf6c87b7e7daf, 0x6045c6bdda77a4d0, 0x4066f541811e1e
};


void print_bytes(unsigned char *p, size_t s)
{
    int i;
    for (i = 0; i < s; ++i) {
        printf("%02x", p[i]);
    }
    printf("\n");
}


unsigned int projective_Z_isAtInfinity(const digit_t* x)
{ 
    unsigned int i;
    
    for (i = 0; i < NWORDS_FIELD; i++) 
    {
        if (x[i] != projective_Z_AtInfinity503[i]) 
        {
            return false;
        }
    }
    return true;
}


void generate_2n_torsionPoints(const f2elm_t *xPA, const f2elm_t *xQA, const f2elm_t *xRA,
                               point_proj_t *torsion_1, point_proj_t *torsion_2,
                               const f2elm_t *A24plus, const f2elm_t *C24, const f2elm_t *A_f2element)
{
    unsigned char *random_k = malloc(NWORDS64_ORDER);
    random_mod_order_A(random_k); //Integer in [0, 2^n -1]

    int baseFound = 0;
    while (baseFound == 0)
    {
        while (1)
        {
            random_mod_order_A(random_k);       //Integer in [0, 2^n -1] - More entropy
            mp_shiftr1(random_k, NWORDS_ORDER); //[0, 2^(n-1) -1]
            mp_shiftl1(random_k, NWORDS_ORDER); //[0, 2^n -2] with only EVEN numbers by dropping the LSB

            LADDER3PT(xPA, xQA, xRA, (digit_t *)random_k, ALICE, torsion_1, A_f2element); //torsion_1 = P + random_k1*Q
            //Now torsion_1 is guaranteed to have order 2^n

            //Same as above
            random_mod_order_A(random_k);
            mp_shiftr1(random_k, NWORDS_ORDER); //[0, 2^(n-1) -1]
            mp_shiftl1(random_k, NWORDS_ORDER); //[0, 2^n -2] with only EVEN numbers by dropping the LSB

            f2elm_t xQminusBase1;
            fp2sub(xQA, torsion_1[0], xQminusBase1);
            LADDER3PT(xPA, xQA, xRA, (digit_t *)random_k, ALICE, torsion_2, A_f2element); //torsion_2 = P + random_k2*Q
            //torsion_2 = torsion_1 + random_k2*Q = (P + random_k1*Q) + random_k2*Q = P + (random_k1 + random_k2)*Q
            //Now torsion_2 is guaranteed to have order 2^n

            /* - Begin check with Tate/Weil Pairing - */

            point_proj_t copy_T1;
            fp2copy(xPA[0], copy_T1->X);
            fp2copy(xPA[1], copy_T1->Z);

            point_proj_t copy_T2;
            fp2copy(xQA[0], copy_T2->X);
            fp2copy(xQA[1], copy_T2->Z);

            xDBLe(torsion_1, copy_T1, A24plus, C24, OALICE_BITS - 2);
            xDBLe(torsion_2, copy_T2, A24plus, C24, OALICE_BITS - 2);

            //Check (affine) x0 != x1 via projective Cross-Multiplication X0*Z1 != X1*Z0
            f2elm_t xT1_zT2, xT2_zT1;
            fp2mul_mont(&copy_T1[0], &copy_T2[1], xT1_zT2);
            fp2mul_mont(&copy_T1[1], &copy_T2[0], xT2_zT1);

            if ((xT1_zT2[0] == xT2_zT1[0]) && (xT1_zT2[1] == xT2_zT1[1]))
            {
                printf("Error: linear dependent points for 2^n basis!\n");
                break;
            }

            //Check torsion_1 order 2^(n-1)
            xDBL(copy_T1, copy_T1, A24plus, C24);

            if (projective_Z_isAtInfinity(copy_T1->Z))
            {
                printf("Error: [2^(n-1)]torsion_1 is the Point at Infinity!\n");
                break;
            }

            //Check torsion_1 order 2^(n)
            xDBL(copy_T1, copy_T1, A24plus, C24);

            if (!projective_Z_isAtInfinity(copy_T1->Z))
            {
                printf("Error: [2^n]torsion_1 is NOT the Point at Infinity!\n");
                break;
            }

            //Check torsion_2 order 2^(n-1)
            xDBL(copy_T2, copy_T2, A24plus, C24);

            if (projective_Z_isAtInfinity(copy_T2->Z))
            {
                printf("Error: [2^(n-1)]torsion_2 is the Point at Infinity!\n");
                break;
            }

            //Check torsion_2 order 2^(n)
            xDBL(copy_T2, copy_T2, A24plus, C24);

            if (!projective_Z_isAtInfinity(copy_T2->Z))
            {
                printf("Error: [2^n]torsion_2 is NOT the Point at Infinity!\n");
                break;
            }

            //Everything is OK
            baseFound = 1;
            break;
        }
    }

    free(random_k); //Free memory for the random integer
}

//Computes the 2n isogenous curve 
//Inputs = isogenyKernel,
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_2n_isogenousCurve_Only(const point_proj_t* isogenyKernel,
                                const f2elm_t* A24plus, const f2elm_t* C24,
                                f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
)
{
    point_proj_t pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeffForPointIsogenies[3];
    unsigned int i, row, m, index, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

        xDBLe(isogenyKernel, S, isogen_A24plus, isogen_C24, (int)(OALICE_BITS-1));
        get_2_isog(S, isogen_A24plus, isogen_C24);
        eval_2_isog(isogenyKernel, S);
#endif

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(isogenyKernel, isogenyKernel, isogen_A24plus, isogen_C24, (int)(2*m));
            index += m;
        }
        get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);

        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        eval_4_isog(isogenyKernel, coeffForPointIsogenies);

        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
}


//Computes the 2n isogenous curve along with the triple phi(P), phi(Q), phi(R)
//Inputs = isogenyKernel,
//         secondBasePoint (point to "isogenize")
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = phiP, phiQ, phiR, phi_BaseT (final points)
//          isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_2n_isogenousCurve_4pt_CurveCoeff(const point_proj_t* isogenyKernel, const point_proj_t* secondBasePoint,
                           point_proj_t* phiP, point_proj_t* phiQ, point_proj_t* phiR, point_proj_t* phi_BaseT, //Outputs
                           const f2elm_t* A24plus, const f2elm_t* C24,
                           f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
)
{
    point_proj_t pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeffForPointIsogenies[3];
    unsigned int i, row, m, index, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    init_basis((digit_t*)B_gen, phiP[0], phiQ[0], phiR[0]);
    fpcopy((digit_t*)&Montgomery_one, phiP[1]->Z);
    fpcopy((digit_t*)&Montgomery_one, phiQ[1]->Z);
    fpcopy((digit_t*)&Montgomery_one, phiR[1]->Z);

    fp2copy(secondBasePoint[0], phi_BaseT[0]);
    fp2copy(secondBasePoint[1], phi_BaseT[1]);

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

        xDBLe(isogenyKernel, S, isogen_A24plus, isogen_C24, (int)(OALICE_BITS-1));
        get_2_isog(S, isogen_A24plus, isogen_C24);
        eval_2_isog(phiP, S);
        eval_2_isog(phiQ, S);
        eval_2_isog(phiR, S);
        eval_2_isog(phi_BaseT, S);
        eval_2_isog(isogenyKernel, S);
#endif

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(isogenyKernel, isogenyKernel, isogen_A24plus, isogen_C24, (int)(2*m));
            index += m;
        }
        get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);

        for (i = 0; i < npts; i++)
        {
            eval_4_isog(pts[i], coeffForPointIsogenies);
        }
        eval_4_isog(phiP, coeffForPointIsogenies);
        eval_4_isog(phiQ, coeffForPointIsogenies);
        eval_4_isog(phiR, coeffForPointIsogenies);
        eval_4_isog(phi_BaseT, coeffForPointIsogenies);

        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    eval_4_isog(phiP, coeffForPointIsogenies);
    eval_4_isog(phiQ, coeffForPointIsogenies);
    eval_4_isog(phiR, coeffForPointIsogenies);
    eval_4_isog(phi_BaseT, coeffForPointIsogenies);

    inv_3_way(phiP[1], phiQ[1], phiR[1]);
    fp2inv_mont(phi_BaseT[1]);
    fp2mul_mont(phiP[0], phiP[1], phiP[0]);
    fp2mul_mont(phiQ[0], phiQ[1], phiQ[0]);
    fp2mul_mont(phiR[0], phiR[1], phiR[0]);
    fp2mul_mont(phi_BaseT[0], phi_BaseT[1], phi_BaseT[0]);
}


//Computes the 2n isogenous curve along with the triple phi(P), phi(Q), phi(R)
//Inputs = isogenyKernel,
//         P, Q, R (initial given points to "isogenize")
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = phiP, phiQ, phiR (final points)
void get_2n_isogenousCurve_3pt(const point_proj_t* isogenyKernel,
                                const point_proj_t* P, const point_proj_t* Q, const point_proj_t* R, //Inputs
                                point_proj_t* phiP, point_proj_t* phiQ, point_proj_t* phiR, //Outputs
                                const f2elm_t* A24plus, const f2elm_t* C24 //Inputs
                                )
{
    point_proj_t pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeffForPointIsogenies[3], isogen_A24plus, isogen_C24;
    unsigned int i, row, m, index, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    fp2copy(P[0], phiP[0]);
    fp2copy(P[1], phiP[1]);
    fp2copy(Q[0], phiQ[0]);
    fp2copy(Q[1], phiQ[1]);
    fp2copy(R[0], phiR[0]);
    fp2copy(R[1], phiR[1]);

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

        xDBLe(isogenyKernel, S, isogen_A24plus, isogen_C24, (int)(OALICE_BITS-1));
        get_2_isog(S, isogen_A24plus, isogen_C24);
        eval_2_isog(phiP, S);
        eval_2_isog(phiQ, S);
        eval_2_isog(phiR, S);
        eval_2_isog(isogenyKernel, S);
#endif

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(isogenyKernel, isogenyKernel, isogen_A24plus, isogen_C24, (int)(2*m));
            index += m;
        }
        get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);

        eval_4_isog(phiP, coeffForPointIsogenies);
        eval_4_isog(phiQ, coeffForPointIsogenies);
        eval_4_isog(phiR, coeffForPointIsogenies);

        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        eval_4_isog(isogenyKernel, coeffForPointIsogenies);

        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    eval_4_isog(phiP, coeffForPointIsogenies);
    eval_4_isog(phiQ, coeffForPointIsogenies);
    eval_4_isog(phiR, coeffForPointIsogenies);

    inv_3_way(phiP[1], phiQ[1], phiR[1]);
    fp2mul_mont(phiP[0], phiP[1], phiP[0]);
    fp2mul_mont(phiQ[0], phiQ[1], phiQ[0]);
    fp2mul_mont(phiR[0], phiR[1], phiR[0]);
}


//Computes the 3m isogenous curve only.
//This function is the first isogeny computed by Bob according to the SIDH-based OT protocol
//Inputs = isogenyKernel,
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_3m_isogenousCurve_Only(const point_proj_t* isogenyKernel,
                                const f2elm_t* A24plus, const f2elm_t* C24,
                                f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
                                )
{
    point_proj_t phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeffForPointIsogenies[3];
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    unsigned char ct_uncomp[UNCOMPRESSEDPK_BYTES] = {0};

    // Initialize basis points
    init_basis((digit_t*)B_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);

    // Traverse tree
    //index = 0;
    for (row = 1; row < MAX_Bob; row++)
    {
        while (index < MAX_Bob-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(isogenyKernel, isogenyKernel, isogen_C24, isogen_A24plus, (int)m);
            index += m;
        }
        get_3_isog(isogenyKernel, isogen_C24, isogen_A24plus, coeffForPointIsogenies);
        
        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        eval_3_isog(isogenyKernel, coeffForPointIsogenies); 
        
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_3_isog(isogenyKernel, isogen_C24, isogen_A24plus, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
}


//Computes the coefficients (x, y) as Alice_phi(T_i) = x*U_i + y*V_i
void linearCombination_Random(const int Alice_Secrets,
                              const point_proj_t* Bob_2n_Basis, const point_proj_t* Alice_Base_phiT,
                              const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                              unsigned char* x_y
                              )
{
    point_proj_t xU, difference, xU_plus_yV;
    unsigned char* random_x = malloc(NWORDS64_ORDER);
    unsigned char* random_y = malloc(NWORDS64_ORDER);

    for(int i = 0; i < Alice_Secrets; i++)
    {
        //Put these inside the do{}while()
        //From here
        random_mod_order_A(random_x);               //x = random
        random_mod_order_A(random_y);               //y = random
        x_y[2*i*SECRETKEY_A_BYTES] = random_x;      //Assign it to x_y array
        x_y[(2*i*SECRETKEY_A_BYTES) +1] = random_y; //Assign it to x_y array

        pointMultiplication(random_x, Bob_2n_Basis[3*i], curveCoeff_A, A24plus, C24, xU);    //xU = random_x * Bob_2n_Basis[3*i]
        pointDifference(Bob_2n_Basis[3*i +1], xU, curveCoeff_A, A24plus, C24, difference);  //difference = Bob_2n_Basis[3*i +1] - xU = V - xU (paper style)
        LADDER3PT(xU->X, Bob_2n_Basis[3*i +1]->X, difference->X, random_y, ALICE, xU_plus_yV, curveCoeff_A); //xU_plus_yV = x*U + y*V (paper style)
        //to here

        //do
        //{
            //Put those lines above, here
        //}
        //while((xU_plus_yV->X != Alice_Base_phiT[2*i]->X) && (xU_plus_yV->Z != Alice_Base_phiT[2*i]->Z));            //while(xU_plus_yV) != phi(T_i)
    }
}


//Computes the multiplication k*P defined in Montgomery projective notation
void pointMultiplication(const digit_t* k, const point_proj_t* P,
                         const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                         point_proj_t* kP)
{
    f2elm_t zero = {0, 0};
    
    fp2copy(zero, kP[0]);
    fp2copy(zero, kP[1]); //kP = Point at infinity

    for(int i = 8*sizeof(k); i >= 0; i--)
    {
        if(((int)k >> i) & 0x01) //Take the i-th bit
        {
            pointAddition(kP, P, curveCoeff_A, A24plus, C24, kP); //kP = kP + B
            xDBL(P, P, A24plus, C24);                             //B = 2B
        }
        else
        {
            pointAddition(kP, P, curveCoeff_A, A24plus, C24, P);  //B = kP + B
            xDBL(kP, kP, A24plus, C24);                           //kP = 2A
        }
    }
}


//Sums two points defined in Montgomery projective notation
void pointAddition(const point_proj* P, const point_proj* Q,
                   const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                   point_proj* sum)
{
    if((P->X == Q->X) && (P->Z == Q->Z)) //if(A == B)
    {
        fp2copy(Q->X, sum->X);
        fp2copy(Q->Z, sum->Z);
        xDBL(P, sum, A24plus, C24);
        return;
    }

    if(is_felm_zero(P->Z[0]) && is_felm_zero(P->Z[1])) //P is the Point at Infinity
    {
        fp2copy(Q->X, sum->X);
        fp2copy(Q->Z, sum->Z);
        return;
    }

    if(is_felm_zero(Q->Z[0]) && is_felm_zero(Q->Z[1])) //Q is the Point at Infinity
    {
        fp2copy(P->X, sum->X);
        fp2copy(P->Z, sum->Z);
        return;
    }

    projectiveSum(P, Q, curveCoeff_A, sum);
}


//Computes C = Q - P in projective coordinates
void pointDifference(const point_proj* Q, const point_proj* P,
                     const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                     point_proj* difference)
{
    point_proj P_neg;
    fp2copy(P->X, P_neg.X);
    fp2copy(P->Z, P_neg.Z);
    fpneg(P_neg.Z); //Negate the Z only

    pointAddition(Q, &P_neg, curveCoeff_A, A24plus, C24, difference);
}


//Sums two points using 15M, 6A
//Formule: https://www.nayuki.io/page/elliptic-curve-point-addition-in-projective-coordinates
void projectiveSum(const point_proj* A, const point_proj* B, const f2elm_t* curveCoeff_A, point_full_proj* pointSUM)
{
    f2elm_t yA, yB;
    getProjective_Y(A, curveCoeff_A, yA);
    getProjective_Y(B, curveCoeff_A, yB);

    f2elm_t t0, t1, T, T_squared, u0, u1, U, U_squared, U_cube, V, W, final_X, final_Y, final_Z;

    fp2mul_mont(yA, B->Z, t0);  //t0 = yA*zB
    fp2mul_mont(yB, A->Z, t1);  //t1 = yB*zA
    fp2sub(t0, t1, T);          //T = t0-t1

    fp2mul_mont(A->X, B->Z, u0);  //u0 = xA*zB
    fp2mul_mont(B->X, A->Z, u1);  //u1 = xB*zA
    fp2sub(u0, u1, U);            //U = u0-u1

    fp2mul_mont(A->Z, B->Z, V);  //V = zA*zB

    fp2sqr_mont(T, T_squared);
    fp2sqr_mont(U, U_squared);
    fp2mul_mont(U, U_squared, U_cube);

    fp2mul_mont(T_squared, V, W);
    fp2add(u0, u1, t1);
    fp2mul_mont(t1, U_squared, t1);
    fp2sub(W, t1, W); //W = T^2*V - U^2(u0+u1)

    fp2mul_mont(U, W, final_X);         //X = U*W
    fp2mul_mont(U_cube, V, final_Z);    //Z = U^3*V

    fp2mul_mont(t0, U_cube, t0);
    fp2mul_mont(u0, U_squared, u0);
    fp2sub(u0, W, final_Y);
    fp2mul_mont(final_Y, T, final_Y);
    fp2sub(final_Y, t0, final_Y);       //Y = T(u0*U^2 - W) - t0*U^3

    fp2copy(final_X, pointSUM->X);
    fp2copy(final_Y, pointSUM->Y);
    fp2copy(final_Z, pointSUM->Z);
}


//Computes the projective Y coordinate of a given projective point P=(X:Z)
void getProjective_Y(const point_proj* P, const f2elm_t* curveCoeff_A, f2elm_t* projective_Y)
{
    f2elm_t temp0;
    fp2sqr_mont(P->X, temp0);                        //temp0 = x^2
    fp2mul_mont(temp0, curveCoeff_A, projective_Y);  //projective_Y = curveCoeff_A * x^2
    fp2mul_mont(P->X, temp0, temp0);                 //temp0 = x^3
    fp2add(temp0, projective_Y, temp0);              //temp0 = x^3 + curveCoeff_A * x^2
    fp2add(temp0, P->X, temp0);                      //temp0 = x^3 + curveCoeff_A * x^2 + x
    sqrt_Fp2(temp0, projective_Y);                   //projective_Y = sqrt(x^3 + curveCoeff_A * x^2 + x)
}


//Convert an Affine point to its Projective form assuming the affine point is normalised
void convertAffine2Projective(point_t* affine, point_proj_t* projective)
{
    f2elm_t one = {0};
    fpcopy((digit_t*)&Montgomery_one, one[0]);

    fp2copy(affine[0], projective[0]);
    fp2copy(one, projective[1]);
}


//Convert an Affine point to its Projective form assuming the affine point is normalised
void convertAffine2_FullProjective(point_t* affine, point_full_proj_t* projective)
{
    f2elm_t one = {0};
    fpcopy((digit_t*)&Montgomery_one, one[0]);

    fp2copy(affine[0], projective[0]->X);
    fp2copy(affine[1], projective[0]->Y);
    fp2copy(one, projective[0]->Z);
}


//Generate a 2^n torsion entangled basis in projective coordinates
void generate_2n_entangledBasis(const f2elm_t* A, point_proj_t* B1, point_proj_t* B2)
{
    point_t S1, S2;
    unsigned char bit, r = 0;
    get_2_torsion_entangled_basis_compression(A, S1, S2, &bit, &r);

    convertAffine2Projective(S1, B1);
    convertAffine2Projective(S1, B1);
}


//Generate a 2^n torsion entangled basis in projective coordinates
void generate_2n_FULL_entangledBasis(const f2elm_t* A, point_full_proj_t* B1, point_full_proj_t* B2)
{
    point_t S1, S2;
    unsigned char bit, r = 0;
    get_2_torsion_entangled_basis_compression(A, S1, S2, &bit, &r);

    convertAffine2_FullProjective(S1, B1);
    convertAffine2_FullProjective(S1, B1);
}


//Converts a (X:Y:Z) point to (X:Z) form
void Convert_Full_to_Short(point_full_proj_t* input, point_proj_t* output)
{
    fp2copy(input[0], output[0]);
    fp2copy(input[2], output[1]);
}


//Computes the Weil Pairing using Tate algorithm
void compute_Tate_Order2n(const point_full_proj_t* P, const point_full_proj_t* Q, f2elm_t* a, const f2elm_t* A_f2element, int t, f2elm_t* n)
{
    point_proj_t P_XZ, Q_XZ, P_XZ_2n, Q_XZ_2n;
    point_full_proj_t P_XYZ, Q_XYZ;
    f2elm_t one_fp2 = {0}, two = {0}, A2, A24, tmp_Y;

    Convert_Full_to_Short(P, P_XZ);
    Convert_Full_to_Short(Q, Q_XZ);

    fpcopy((digit_t*)&Montgomery_one, one_fp2[0]);
    fp2div2(A_f2element,A2);
    fpcopy(one_fp2[0],two[0]);
    fpadd(two[0], two[0], two[0]);
    fp2add(A_f2element,two,A24);
    fp2div2(A24,A24);
    fp2div2(A24,A24);

    Double(P_XZ, P_XZ_2n, A24, OALICE_BITS-1); //P_XZ_2n = [2^(n-1)]P_XZ
    Double(Q_XZ, Q_XZ_2n, A24, OALICE_BITS-1); //Q_XZ_2n = [2^(n-1)]Q_XZ

    getProjective_Y(P_XZ_2n, A_f2element, tmp_Y);
    fp2copy(P_XZ_2n[0].X, P_XYZ[0].X);
    fp2copy(tmp_Y,        P_XYZ[0].Y);
    fp2copy(P_XZ_2n[0].Z, P_XYZ[0].Z);

    getProjective_Y(Q_XZ_2n, A_f2element, tmp_Y);
    fp2copy(Q_XZ_2n[0].X, Q_XYZ[0].X);
    fp2copy(tmp_Y,        Q_XYZ[0].Y);
    fp2copy(Q_XZ_2n[0].Z, Q_XYZ[0].Z);

    Tate_pairings_2_torsion(P_XYZ, Q_XYZ, a, t, n); //With P_XZ_2n and Q_XZ_2n we only need an order 2 pairing instead of an order 2^n
}


//Encode and Decode a secret using a j_invariant as seed for a KDF
void Enc(unsigned char* secret, f2elm_t j_invariant, unsigned char *XORed_secret)
{ 
    unsigned char charEncoded_j_inv[FP2_ENCODED_BYTES] = {0};
    unsigned char seed[MSG_BYTES];

    fp2_encode(j_invariant, charEncoded_j_inv);
    shake256(seed, MSG_BYTES, charEncoded_j_inv, FP2_ENCODED_BYTES); //seed = KFD(j_inv)

    for (int i = 0; i < MSG_BYTES; i++) 
    {
        XORed_secret[i] = secret[i] ^ seed[i];
    }
}

int64_t cpucycles()
{ 
    unsigned int hi, lo;

    asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
}


//SIOT




//Computes the 2n isogenous curve along with the triple phi(P), phi(Q), phi(R)
//Inputs = isogenyKernel,
//         P, Q, R (initial given points to "isogenize")
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = phiP, phiQ, phiR (final points)
//          isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_2n_isogenousCurve_3pt_CurveCoeff(const point_proj_t* isogenyKernel,
                                        const point_proj_t* P, const point_proj_t* Q, const point_proj_t* R, //Inputs
                                        point_proj_t* phiP, point_proj_t* phiQ, point_proj_t* phiR, //Outputs
                                        const f2elm_t* A24plus, const f2elm_t* C24, //Inputs
                                        f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
)
{
    point_proj_t pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeffForPointIsogenies[3];
    unsigned int i, row, m, index, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    fp2copy(P[0], phiP[0]);
    fp2copy(P[1], phiP[1]);
    fp2copy(Q[0], phiQ[0]);
    fp2copy(Q[1], phiQ[1]);
    fp2copy(R[0], phiR[0]);
    fp2copy(R[1], phiR[1]);

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);
    
    
#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

        xDBLe(isogenyKernel, S, isogen_A24plus, isogen_C24, (int)(OALICE_BITS-1));
        get_2_isog(S, isogen_A24plus, isogen_C24);
        eval_2_isog(phiP, S);
        eval_2_isog(phiQ, S);
        eval_2_isog(phiR, S);
        eval_2_isog(isogenyKernel, S);
#endif

    // Traverse tree
    index = 0;
    for (row = 1; row < MAX_Alice; row++)
    {
        while (index < MAX_Alice-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(isogenyKernel, isogenyKernel, isogen_A24plus, isogen_C24, (int)(2*m));
            index += m;
        }
        get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);

        for (i = 0; i < npts; i++)
        {
            eval_4_isog(pts[i], coeffForPointIsogenies);
        }
        eval_4_isog(phiP, coeffForPointIsogenies);
        eval_4_isog(phiQ, coeffForPointIsogenies);
        eval_4_isog(phiR, coeffForPointIsogenies);

        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(isogenyKernel, isogen_A24plus, isogen_C24, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    eval_4_isog(phiP, coeffForPointIsogenies);
    eval_4_isog(phiQ, coeffForPointIsogenies);
    eval_4_isog(phiR, coeffForPointIsogenies);

    inv_3_way(phiP[1], phiQ[1], phiR[1]);
    fp2mul_mont(phiP[0], phiP[1], phiP[0]);
    fp2mul_mont(phiQ[0], phiQ[1], phiQ[0]);
    fp2mul_mont(phiR[0], phiR[1], phiR[0]);
}



//Computes the 3m isogenous curve.
//Inputs = isogenyKernel,
//         P, Q, R (initial given points to "isogenize")
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = phiP, phiQ, phiR (final points)
//          isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_3m_isogenousCurve_3pt_Coeff(const point_proj_t* isogenyKernel,
                                    const point_proj_t* P, const point_proj_t* Q, const point_proj_t* R, //Inputs
                                    const f2elm_t* A24plus, const f2elm_t* C24,
                                    point_proj_t* phiP, point_proj_t* phiQ, point_proj_t* phiR, //Outputs
                                    f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
                                )
{
    point_proj_t pts[MAX_INT_POINTS_BOB];
    f2elm_t coeffForPointIsogenies[3];
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    unsigned char ct_uncomp[UNCOMPRESSEDPK_BYTES] = {0};

    fp2copy(P[0], phiP[0]);
    fp2copy(P[1], phiP[1]);
    fp2copy(Q[0], phiQ[0]);
    fp2copy(Q[1], phiQ[1]);
    fp2copy(R[0], phiR[0]);
    fp2copy(R[1], phiR[1]);

    fp2copy(A24plus, isogen_A24plus);
    fp2copy(C24, isogen_C24);

    // Traverse tree
    //index = 0;
    for (row = 1; row < MAX_Bob; row++)
    {
        while (index < MAX_Bob-row)
        {
            fp2copy(isogenyKernel[0], pts[npts]->X);
            fp2copy(isogenyKernel[1], pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(isogenyKernel, isogenyKernel, isogen_C24, isogen_A24plus, (int)m);
            index += m;
        }
        get_3_isog(isogenyKernel, isogen_C24, isogen_A24plus, coeffForPointIsogenies);

        eval_3_isog(phiP, coeffForPointIsogenies);
        eval_3_isog(phiQ, coeffForPointIsogenies);
        eval_3_isog(phiR, coeffForPointIsogenies);
        
        fp2copy(pts[npts-1]->X, isogenyKernel[0]);
        fp2copy(pts[npts-1]->Z, isogenyKernel[1]);
        eval_3_isog(isogenyKernel, coeffForPointIsogenies); 
        
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_3_isog(isogenyKernel, isogen_C24, isogen_A24plus, coeffForPointIsogenies);
    fp2add(isogen_A24plus, isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    fp2div2(isogen_C24, isogen_C24);
    eval_3_isog(phiP, coeffForPointIsogenies);
    eval_3_isog(phiQ, coeffForPointIsogenies);
    eval_3_isog(phiR, coeffForPointIsogenies);

    inv_3_way(phiP[1], phiQ[1], phiR[1]);
    fp2mul_mont(phiP[0], phiP[1], phiP[0]);
    fp2mul_mont(phiQ[0], phiQ[1], phiQ[0]);
    fp2mul_mont(phiR[0], phiR[1], phiR[0]);
}

//Computes the curve constant C from A24plus and C24 constants
//Source: SIKE - Supersingular Isogeny Key Encapsulation, https://sike.org/ , April 17, 2019
//Appendix A - Explicit algorithms for isogen‘ and isoex‘: Optimized implementation
// (A+ 24 : C24) ∼ (A + 2C : 4C)
void get_A2_from_A24p_C24(const f2elm_t* A24plus, const f2elm_t* C24, f2elm_t* A_f2element)
{
    fp2add(C24, C24, A_f2element); // A = 2C
    fp2sub(A_f2element, A24plus, A_f2element); // A = A - A+24 = 2C - A+24
}


//Check whether the decrypted message is equal to the initial secret
int check_Result(unsigned char *secret, unsigned char *decrypted, size_t message_length)
{
    int i, check;
    for (i = 0; i < message_length; ++i) 
    {
        if(secret[i] != decrypted[i])
        {
            return check = 0;
        }
    }
    return check = 1;
}