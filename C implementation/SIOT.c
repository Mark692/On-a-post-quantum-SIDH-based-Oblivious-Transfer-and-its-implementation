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

#include "obliviousTransfer.c"
#include "SIOT.h"



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
    eval_3_isog(phiP, coeffForPointIsogenies);
    eval_3_isog(phiQ, coeffForPointIsogenies);
    eval_3_isog(phiR, coeffForPointIsogenies);

    inv_3_way(phiP[1], phiQ[1], phiR[1]);
    fp2mul_mont(phiP[0], phiP[1], phiP[0]);
    fp2mul_mont(phiQ[0], phiQ[1], phiQ[0]);
    fp2mul_mont(phiR[0], phiR[1], phiR[0]);
}
