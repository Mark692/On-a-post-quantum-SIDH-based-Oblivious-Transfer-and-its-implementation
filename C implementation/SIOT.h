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

   

#ifndef SIOT_H
#define SIOT_H



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
);


//Computes the 3m isogenous curve only.
//This function is the first isogeny computed by Bob according to the SIDH-based OT protocol
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
);

#endif //SIOT_H
