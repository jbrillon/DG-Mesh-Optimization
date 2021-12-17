#ifndef VAR_H
#define VAR_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>

/* Setup parameters */
extern std::string problem, subproblem, namesetup;
extern int P, Nel, np, Mdim, npG;
extern double T0, Tf, CFL, TIME_STEP, DELTA_X_MIN;

/* Vectors */
extern double *z, *w, *JacobPolyStore, *elBounds, *xDomainBounds, *Jac, *phi1,
 *phi2, *dphi1, *dphi2, *xG, *D_flat, *Dt_flat;

/* Matrices */
extern double **DUMMY_INNER_PROD, **DUMMY_MAT, **u_exact, **u_delta, **u_hat, **u0, **x, **f, **f_hat, **f_star, **D, **Dt;
extern int **map, **map_frequency, **map_physical;
extern double **phi_L, **phi_R, **RHS_VEC, **u_dummy, **q_source, **q_phi_innerproduct;

/* Blocks */
extern double ***M, ***C, ***L;

/* LAPACK VAR */
extern unsigned char TRANS;
extern int NRHS, INFO, *ipiv, *ipiv_global;
extern double *DUMMY_VEC, *DUMMY_MAT_FLAT;
extern double **RK_INPUT_VEC, **K1, **K2, **K3, **K4;

/* Physical and mathematical constants */
extern double PI, LinAdvSpeed;

/* Transient */
extern int store_transient, iout, DUMP_INTERVAL;

/* Error estimates */
extern int nDOF, ConvStudyFlag, ConvStudyType;
extern double *epsG;

/* MDO */
extern int arctan_NumShocks;
extern double arctan_const, arctan_magnitude;
extern double *arctan_ShockMag, *arctan_ShockLoc;

/* Implicit */
extern int GlobalDim;
extern double *u_hat_global, *RHS_VEC_global, *q_phi_innerproduct_global, *LHS_MAT_global_FLAT;
extern double **LHS_MAT_global, **M_global, **C_global, **L_global, **PHI_LL_MAT, **PHI_LR_MAT, **PHI_MAT_global;

#endif