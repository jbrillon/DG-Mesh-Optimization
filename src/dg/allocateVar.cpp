#include <iostream>
#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <stdio.h>
// My standard libraries
#include "allocatePointersLib.h"
// My SEM libraries
#include "var.h"
#include "allocateVar.h"

using namespace std;

void allocateVariables()
{
    //-----------------------------------------------
    //          Allocate all
    //-----------------------------------------------
    /* Vectors */
    ipiv = ivector(Mdim); // for LAPACK
    DUMMY_VEC = dvector(Mdim); // for LAPACK
    DUMMY_MAT_FLAT = dvector(Mdim*Mdim); // for LAPACK
    z = dvector(np); // zeros (xi) (standard domain)
    w = dvector(np); // weights
    JacobPolyStore = dvector(np);
    elBounds = dvector(Nel+1);
    Jac = dvector(Nel);
    phi1 = dvector(np);
    phi2 = dvector(np);
    dphi1 = dvector(np);
    dphi2 = dvector(np);
    xG = dvector(npG);
    D_flat = dvector(np*np);
    Dt_flat = dvector(np*np);
    /* Matrices */
    phi_L = dmatrix(Nel,Mdim); // should change to dvector
    phi_R = dmatrix(Nel,Mdim); // should change to dvector
    RHS_VEC = dmatrix(Nel,Mdim);
    DUMMY_INNER_PROD = dmatrix(Nel, Mdim);
    DUMMY_MAT = dmatrix(Mdim, Mdim);
    u_exact = dmatrix(Nel, np); // exact solution in physical domain
    u_delta = dmatrix(Nel, np); // solution in physical domain
    u_dummy = dmatrix(Nel, np); // solution (dummy storage) in physical domain
    q_source = dmatrix(Nel, np); // source term
    q_phi_innerproduct = dmatrix(Nel, np); // inner product of source term and phi
    u_hat = dmatrix(Nel, Mdim); // weights for each element
    RK_INPUT_VEC = dmatrix(Nel, Mdim); // for RK4
    K1 = dmatrix(Nel, Mdim); // for RK4
    K2 = dmatrix(Nel, Mdim); // for RK4
    K3 = dmatrix(Nel, Mdim); // for RK4
    K4 = dmatrix(Nel, Mdim); // for RK4
    u0 = dmatrix(Nel, np); // initial condition in physical space
    x = dmatrix(Nel, np); // element domain stored
    f = dmatrix(Nel, np); // physical flux in physical domain
    f_hat = dmatrix(Nel, Mdim); // flux weights for each element
    f_star = dmatrix(Nel, 2); // numerical flux in physical domain
    D = dmatrix(np, np); // dgll diff. matrix
    Dt = dmatrix(np, np); // transpose of D
    map = imatrix(Nel, Mdim); // mapping matrix
    map_frequency = imatrix(Nel, Mdim); // mapping matrix: frequency space
    map_physical = imatrix(Nel, np); // mapping matrix: physical space
    /* Blocks */
    M = dblock(Mdim, Mdim, Nel); // mass matrix of element
    C = dblock(Mdim, Mdim, Nel); // d/dx matrix of element
    L = dblock(Mdim, Mdim, Nel); // d2/dx2 matrix of element
    //-----------------------------------------------
    //          Error estimates
    //-----------------------------------------------
    epsG = dvector(npG);
    //-----------------------------------------------
    //          Global / Implicit Stuff
    //-----------------------------------------------
    ipiv_global = ivector(GlobalDim); // for LAPACK
    u_hat_global = dvector(GlobalDim);
    RHS_VEC_global = dvector(GlobalDim);
    q_phi_innerproduct_global = dvector(GlobalDim);
    LHS_MAT_global_FLAT = dvector(GlobalDim*GlobalDim);
    LHS_MAT_global = dmatrix(GlobalDim, GlobalDim);
    M_global = dmatrix(GlobalDim, GlobalDim);
    C_global = dmatrix(GlobalDim, GlobalDim);
    L_global = dmatrix(GlobalDim, GlobalDim);
    PHI_LL_MAT = dmatrix(Mdim, Mdim);
    PHI_LR_MAT = dmatrix(Mdim, Mdim);
    PHI_MAT_global = dmatrix(GlobalDim, GlobalDim);
}