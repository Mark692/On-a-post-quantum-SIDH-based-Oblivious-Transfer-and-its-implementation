//
// Created by markarol on 01/12/19.
//

#include "unUsed.h"

// Convert a Weiestrass (affine) Point into a Montgomery Point on the curve y^2 = x^3 + A*x^2 + x
point_full_proj_t convert_Weier2Mont(const point_full_proj_t pointWeier, const f2elm_t curve_A)
{
    f2elm_t zeroElem = {0}, oneElem = {0}, temp;
    point_full_proj_t montgomeryPoint;

    fpcopy((digit_t*)&Montgomery_one, oneElem[0]);

    if(is_felm_zero(pointWeier->X) && is_felm_zero(pointWeier->Z))
    {
        fp2copy(zeroElem, montgomeryPoint->X);
        fp2copy(oneElem,  montgomeryPoint->Y);
        fp2copy(zeroElem, montgomeryPoint->Z);
        return montgomeryPoint;
    }
    fpmul_mont((digit_t*)threeinv, curve_A[0], temp[0]);
    fpmul_mont((digit_t*)threeinv, curve_A[1], temp[1]);

    fp2sub(temp, pointWeier->X, montgomeryPoint->X); //xMont = xWeier - A/3
    fp2copy(pointWeier->Y, montgomeryPoint->Y);
    fp2copy(oneElem, montgomeryPoint->Z);
    return montgomeryPoint;
}


//Convert a projective point into its affine representation
//Formulae: http://hyperelliptic.org/EFD/g1p/auto-montgom.html
point_t convert_Projective2Affine(const point_proj_t* projectivePoint, const f2elm_t* curve_A)
{
    point_t affineP;
    f2elm_t inverse_zA;

    fp2copy(projectivePoint->Z, inverse_zA);
    fp2inv_mont(inverse_zA); //1/zA
    fp2mul_mont(projectivePoint->X, inverse_zA, affineP->x); //affineP->x = X/Z

    //TODO: Dato che affineP->x = X/Z e vale anche affineP->y = Y/Z
    //Problema: se nella rappresentazione proiettiva abbiamo Y=1 allora
    //la affineP->y sarà sempre 1/Z oppure la devo calcolare tramite la seguente formula?
    f2elm_t temp0, temp1;
    fp2sqr_mont(affineP->x, temp0);        //temp0 = x^2
    fp2mul_mont(temp0, curve_A, temp1);    //temp1 = curve_A * x^2
    fp2mul_mont(affineP->x, temp0, temp0); //temp0 = x^3

    fp2add(temp0, temp1, temp0);           //x^3 + curve_A * x^2
    fp2add(temp0, affineP->x, temp0);      //x^3 + curve_A * x^2 + x
    sqrt_Fp2(temp0, affineP->y);           //y = sqrt(x^3 + curve_A * x^2 + x)

    return affineP;
}


//Computes affine sum between two affine points
point_t affineSUM(const point_t* A, const point_t* B, const f2elm_t* curve_A)
{
    point_t sum;
    f2elm_t t0, t1, t2, t3;

    f2psub(B->y, A->y, t0);     //t0 = yB - yA
    fp2sub(B->x, A->x, t1);     //t1 = xB - xA
    fp2copy(t1, t2);            //t2 = t1
    fp2inv_mont(t2);            //t2 = 1/t2
    fp2mul_mont(t0, t2, t2);    //t2 = (yB - yA)/(xB - xA)
    fp2sqr_mont(t2, t3);        //t3 = [(yB - yA)/(xB - xA)]^2

    fp2sub(t3,     curve_A, sum->x);
    fp2sub(sum->x, A->x,    sum->x);
    fp2sub(sum->x, B->x,    sum->x); //x = t3 - curve_A - xA - xB

    fp2add(A->x, A->x, t0); //t0 = 2*xA
    fp2add(t0, B->x, sum->y);
    fp2add(sum->x, curve_A, sum->y);
    fp2mul_mont(sum->y, t2, sum->y);
    fp2mul_mont(t3, t2, t3); //t3 = [(yB - yA)/(xB - xA)]^3
    fp2sub(sum->y, t3, sum->y);
    fp2sub(sum->y, B->y, sum->y); //y = Formula from the link

    return sum;
}
