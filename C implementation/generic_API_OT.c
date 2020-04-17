//
// Created by markarol on 02/12/19.
//
#include <stdio.h>
#include <stdlib.h>
#include "lib/P503/api.h"
#include "lib/config.h"
#include "lib/P503/P503_internal.h"

void print_bytes(unsigned char *p, size_t s)
{
    int i;
    for (i = 0; i < s; ++i) {
        printf("%02x", p[i]);
    }
    printf("\n");
}


int generate_2n_basis(const f2elm_t xPA, const f2elm_t xQA, const f2elm_t xRA,
                       point_proj* base1, point_proj* base2,
                       const f2elm_t A24plus, const f2elm_t C24, const f2elm_t A_f2element)
{
    unsigned char* random_k;
    random_mod_order_A(random_k); //Integer in [0, 2^n -1]
    mp_shiftr1(random_k, NWORDS_ORDER); //[0, 2^(n-1) -1]
    mp_shiftl1(random_k, NWORDS_ORDER); //[0, 2^n -2] with only EVEN numbers by dropping the LSB

    LADDER3PT(xPA, xQA, xRA, (digit_t*) random_k, ALICE, base1, A_f2element); //base1 = P + random_k1*Q
    //Now base1 is guaranteed to have order 2^n

    //Same as above
    random_mod_order_A(random_k);
    mp_shiftr1(random_k, NWORDS_ORDER); //[0, 2^(n-1) -1]
    mp_shiftl1(random_k, NWORDS_ORDER); //[0, 2^n -2] with only EVEN numbers by dropping the LSB

    f2elm_t xQminusBase1;
    fp2sub(xQA, base1->X, xQminusBase1);
    LADDER3PT(base1->X, xQA, xQminusBase1, (digit_t*) random_k, ALICE, base2, A_f2element); //base2 = base1 + random_k2*Q
    //base2 = base1 + random_k2*Q = (P + random_k1*Q) + random_k2*Q = P + (random_k1 + random_k2)*Q
    //Now base2 is guaranteed to have order 2^n


    //Check TatePairing/WeilPairing
    point_proj_t copyB1;
    fp2copy(base1->X, copyB1->X);
    fp2copy(base1->Z, copyB1->Z);

    point_proj_t copyB2;
    fp2copy(base2->X, copyB2->X);
    fp2copy(base2->Z, copyB2->Z);

    xDBLe(base1, copyB1, A24plus, C24, OALICE_BITS-2);
    xDBLe(base2, copyB2, A24plus, C24, OALICE_BITS-2);


    //Check (affine) x0 != x1 via projective Cross-Multiplication X0*Z1 != X1*Z0
    f2elm_t xB1_zB2, xB2_zB1;
    fp2mul_mont(base1->X, base2->Z, xB1_zB2);
    fp2mul_mont(base1->Z, base2->X, xB2_zB1);
    if(xB1_zB2[0] == xB2_zB1[0]) //Real part
    {
        if(xB1_zB2[1] == xB2_zB1[1]) //Imaginary part
        {
            printf("Error: linear dependent points for 2^n basis!");
            return 0;
        }
    }


    xDBL(copyB1, copyB1, A24plus, C24);
    if(copyB1->Z == 0)
    {
        printf("Error: [2^(n-1)]base1 is the Point at Infinity!");
        return 0;
    }

    xDBL(copyB1, copyB1, A24plus, C24);
    if(copyB1->Z != 0)
    {
        printf("Error: [2^n]base1 is NOT the Point at Infinity!");
        return 0;
    }

    xDBL(copyB2, copyB2, A24plus, C24);
    if(copyB2->Z == 0)
    {
        printf("Error: [2^(n-1)]base2 is the Point at Infinity!");
        return 0;
    }

    xDBL(copyB2, copyB2, A24plus, C24);
    if(copyB2->Z != 0)
    {
        printf("Error: [2^n]base2 is NOT the Point at Infinity!");
        return 0;
    }

    //Everything is OK
    return 1;
}


//Sums two points defined in Montgomery projective notation
void pointAddition(const point_proj* P, const point_proj* Q,
        const f2elm_t curveCoeff_A, const f2elm_t A24plus, const f2elm_t C24,
        point_proj* sum)
{
    if((P->X == Q->X) && (P->Z == Q->Z)) //if(A == B)
    {
        fp2copy(Q->X, sum->X);
        fp2copy(Q->Z, sum->Z);
        xDBL(P, sum, A24plus, C24);
        return;
    }

    if(fp2zero(P->Z[0]) && fp2zero(P->Z[1]))
    {
        fp2copy(Q->X, sum->X);
        fp2copy(Q->Z, sum->Z);
        return;
    }

    if(fp2zero(Q->Z[0]) && fp2zero(Q->Z[1]))
    {
        fp2copy(P->X, sum->X);
        fp2copy(P->Z, sum->Z);
        return;
    }

    projectiveSum(P, Q, curveCoeff_A, sum);
}


//Computes C = Q - P in projective coordinates
void pointDifference(const point_proj* Q, const point_proj* P,
        const f2elm_t curveCoeff_A, const f2elm_t A24plus, const f2elm_t C24,
        point_proj* difference)
{
    point_proj* P_neg;
    fp2copy(P->X, P_neg->X);
    fp2copy(P->Z, P_neg->Z);
    fpneg(P_neg->Z); //Negate the Z only

    pointAddition(Q, P_neg, curveCoeff_A, A24plus, C24, difference);
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


//Computes the projective Y coordinate of a given point P
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