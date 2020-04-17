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

   

#ifndef OBLIVIOUSTRANSFER_H
#define OBLIVIOUSTRANSFER_H

#define __LINUX__   // -> lib/config.h defines OS_TARGET
#define __GNUC__    // -> lib/config.h defines COMPILER COMPILER_GCC
#define _AMD64_     // -> lib/config.h defines TARGET,  digit_t...

#define _GENERIC_   // -> lib/config.h defines the implementation


void print_bytes(unsigned char *p, size_t s);

unsigned int projective_Z_isAtInfinity(const digit_t* x);


void generate_2n_torsionPoints(const f2elm_t *xPA, const f2elm_t *xQA, const f2elm_t *xRA,
                               point_proj_t *torsion_1, point_proj_t *torsion_2,
                               const f2elm_t *A24plus, const f2elm_t *C24, const f2elm_t *A_f2element
);


//Computes the 2n isogenous curve 
//Inputs = isogenyKernel,
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_2n_isogenousCurve_Only(const point_proj_t* isogenyKernel,
                                const f2elm_t* A24plus, const f2elm_t* C24,
                                f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
);



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
);


//Computes the 2n isogenous curve along with the triple phi(P), phi(Q), phi(R)
//Inputs = isogenyKernel,
//         P, Q, R (initial given points to "isogenize")
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = phiP, phiQ, phiR (final points)
//          isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_2n_isogenousCurve_3pt(const point_proj_t* isogenyKernel,
                                     const point_proj_t* P, const point_proj_t* Q, const point_proj_t* R, //Inputs
                                     point_proj_t* phiP, point_proj_t* phiQ, point_proj_t* phiR, //Outputs
                                     const f2elm_t* A24plus, const f2elm_t* C24 //Inputs
);


//Computes the 3m isogenous curve only.
//This function is the first isogeny computed by Bob according to the SIDH-based OT protocol
//Inputs = isogenyKernel,
//         A24plus, C24 (initial curve coefficients)
//
//Outputs = isogen_A24plus, isogen_C24 (isogenous curve coefficients)
void get_3m_isogenousCurve_Only(const point_proj_t* isogenyKernel,
                                const f2elm_t* A24plus, const f2elm_t* C24,
                                f2elm_t* isogen_A24plus, f2elm_t* isogen_C24 //Outputs
);


//Computes the coefficients (x, y) as Alice_phi(T_i) = x*U_i + y*V_i
void linearCombination_Random(const int Alice_Secrets,
                              const point_proj_t* Bob_2n_Basis, const point_proj_t* Alice_Base_phiT,
                              const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                              unsigned char* x_y
);


//Computes the multiplication k*P defined in Montgomery projective notation
void pointMultiplication(const digit_t* k, const point_proj_t* P,
                         const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                         point_proj_t* kP
);


//Sums two points defined in Montgomery projective notation
void pointAddition(const point_proj* P, const point_proj* Q,
                   const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                   point_proj* sum
);


//Computes C = Q - P in projective coordinates
void pointDifference(const point_proj* Q, const point_proj* P,
                     const f2elm_t* curveCoeff_A, const f2elm_t* A24plus, const f2elm_t* C24,
                     point_proj* difference
);


//Sums two points using 15M, 6A
//Formule: https://www.nayuki.io/page/elliptic-curve-point-addition-in-projective-coordinates
void projectiveSum(const point_proj* A, const point_proj* B, const f2elm_t* curveCoeff_A, point_full_proj* pointSUM);


//Computes the projective Y coordinate of a given projective point P=(X:Z)
void getProjective_Y(const point_proj* P, const f2elm_t* curveCoeff_A, f2elm_t* projective_Y);


//Convert an Affine point to its Projective form assuming the affine point is normalised
void convertAffine2Projective(point_t* affine, point_proj_t* projective);


//Convert an Affine point to its Projective form assuming the affine point is normalised
void convertAffine2_FullProjective(point_t* affine, point_full_proj_t* projective);


//Generate a 2^n torsion entangled basis in projective coordinates
void generate_2n_entangledBasis(const f2elm_t* A, point_proj_t* B1, point_proj_t* B2);


//Generate a 2^n torsion entangled basis in projective coordinates
void generate_2n_FULL_entangledBasis(const f2elm_t* A, point_full_proj_t* B1, point_full_proj_t* B2);


//Converts a (X:Y:Z) point to (X:Z) form
void Convert_Full_to_Short(point_full_proj_t* input, point_proj_t* output);


//Computes the Weil Pairing using Tate algorithm
void compute_Tate_Order2n(const point_full_proj_t* P, const point_full_proj_t* Q, f2elm_t* a, const f2elm_t* A_f2element, int t, f2elm_t* n);


//Encode and Decode a secret using a j_invariant as seed for a KDF
void Enc(unsigned char* secret, f2elm_t j_invariant, unsigned char *XORed_secret);


//Used to count cpucycles
//Source: library/tests/test_extras.c
//int64_t cpucycles();








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

//Computes the curve constant C from A24plus and C24 constants
void get_A2_from_A24p_C24(const f2elm_t* A24plus, const f2elm_t* C24, f2elm_t* A_f2element);

//Check whether the decrypted message is equal to the initial secret
int check_Result(unsigned char *secret, unsigned char *decrypted, size_t message_length);


#endif //MARCO_OBLIVIOUSTRANSFER_H
