#include <iostream>
#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <stdio.h>
#include "allocatePointersLib.h"
#include "var.h"


/* Variables read in by readFilesDG */
std::string problem, subproblem, namesetup;
int P, Nel; // polynomial order, number of elements
double T0, Tf, CFL;

/* Time advancement variables */
double TIME_STEP, DELTA_X_MIN;

/* Intergers */
int np, Mdim, npG;

/* Vectors */
double *z, *w, *JacobPolyStore, *elBounds, *Jac, *phi1,
 *phi2, *dphi1, *dphi2, *xG, *D_flat, *Dt_flat;
double *xDomainBounds = dvector(2); // allocate here since it's input before allocating all variables

/* Matrices */
double **DUMMY_INNER_PROD, **DUMMY_MAT, **u_exact, **u_delta, **u_hat, **u0, **x, **f, **f_hat, **f_star, **D, **Dt;
int **map, **map_frequency, **map_physical;
double **phi_L, **phi_R, **RHS_VEC, **u_dummy, **q_source, **q_phi_innerproduct;
double **RK_INPUT_VEC, **K1, **K2, **K3, **K4;

/* Blocks */
double ***M, ***C, ***L;

/* LAPACK VAR */
unsigned char TRANS = 'T';
int NRHS=1, INFO=0, *ipiv, *ipiv_global;
double *DUMMY_VEC, *DUMMY_MAT_FLAT;

/* Physical and mathematical constants */
double PI = 3.14159265359;
double LinAdvSpeed;

/* Transient */
int store_transient=-1, iout, DUMP_INTERVAL=10;

/* Error estimates */
int nDOF, ConvStudyFlag=-1, ConvStudyType=-1;
double *epsG;

/* MDO */
int arctan_NumShocks;
double arctan_const, arctan_magnitude;
double *arctan_ShockMag, *arctan_ShockLoc;

/* Implicit */
int GlobalDim;
double *u_hat_global, *RHS_VEC_global, *q_phi_innerproduct_global, *LHS_MAT_global_FLAT;
double **LHS_MAT_global, **M_global, **C_global, **L_global, **PHI_LL_MAT, **PHI_LR_MAT, **PHI_MAT_global;
