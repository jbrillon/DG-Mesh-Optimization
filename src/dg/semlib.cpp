#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// polynomial library from NEKTAR
#include "polylib.h"
// my libraries
#include "flattenLib.h"
#include "var.h"
#include "semlib.h"

using namespace polylib;
using namespace std;

//===============================================
//          LAPACK FUNCTIONS
//===============================================
extern "C" {extern void dgesv_(int *, int *, double (*), int *, int [], double (*), int *, int*);}

//===============================================
/* All code here has been double checked since 
 * revisiting this project */
//===============================================

//===============================================
/* Element bounds */
//-----------------------------------------------
void elementBounds(int Nel, double *xDomainBounds, double *elBounds)
{
   // for MDO project (vertex points)
   for(int i=0; i<(Nel+1); i++)
   {
      elBounds[i] = xDomainBounds[0] + double(i)*(xDomainBounds[1]-xDomainBounds[0])/double(Nel);
   }
}
//===============================================
/* Mapping function and jacobian evaluation */
//-----------------------------------------------
void chi(int np, int eln, double *z, double *elBounds,
 double **x, double *Jac)
{
   // Maps local to global
   // elementBounds = [lower, upper] bound
   // Mapping description: Linear
   for(int i=0; i<np; i++)
   {
      x[eln][i] = 0.5*((1.0-z[i])*elBounds[eln] + (1.0+z[i])*elBounds[eln+1]);
   }
   // Jacobian (d\chi/d\xi i.e. dx/dz)
   Jac[eln] = 0.5*(elBounds[eln+1] - elBounds[eln]);
}
//===============================================
/* Modal p-type expansion basis */
//-----------------------------------------------
void basis(int np, int P, int i, double *z, double *phi)
{
   if(i==0) 
   {
      for(int k=0; k<np; k++)
      {
         phi[k] = 0.5*(1.0-z[k]);
      }
   }
   else if(i==P)
   {
      for(int k=0; k<np; k++)
      {
         phi[k] = 0.5*(1.0+z[k]);
      }
   }
   else
   {
      jacobfd(np, z, JacobPolyStore, NULL, i-1, 1.0, 1.0); // phi stores evaluation

      for(int k=0; k<np; k++)
      {
         phi[k] = 0.25*(1.0-z[k])*(1.0+z[k])*JacobPolyStore[k];
      }
   }
}
//===============================================
/* Generate modal basis vectors at xi_L and xi_R for
   each element */
//-----------------------------------------------
void basis_vec_xLR(int np, int P, double *z, 
   double *phi1, double **phi_L, double **phi_R)
{
   /* INPUT:
    * - z: standard domain (xi)
    * OUTPUT:
    * - vector of phi_p(z) for p in [1,P]
        for each z=zL,zR (-1,1) of each element
   */
   for(int i=0; i<(P+1); i++)
   {
      basis(np, P, i, z, phi1);
      for(int eln=0; eln<Nel; eln++)
      {
         phi_L[eln][i] = phi1[0];
         phi_R[eln][i] = phi1[np-1];
      }
   }
}
//===============================================
/* Collocation differentiation wrt to x */
//-----------------------------------------------
void diff(int np, int eln, double *phi, double **D, double *Jac, double *dphi)
{
   for(int i=0; i<np; i++)
   {
      dphi[i] = 0.0;
      for(int j=0; j<np; j++)
      {
         dphi[i] += D[i][j]*phi[j];
      }
      dphi[i] /= Jac[eln];
   }
}
//===============================================
/* Gaussian quadrature for phi1*phi2 wrt x */
//-----------------------------------------------
double integr(int np, int eln, double *w, double *Jac, double *phi1, double *phi2)
{
   register double sum = 0.0;

   for(int i=0; i<np; i++)
   {
      sum += w[i]*phi1[i]*phi2[i];
   }
   sum *= Jac[eln];
   return sum;
}
//-----------------------------------------------
double integr_VecMat(int np, int eln, double *w, double *Jac, double *phi1, double **phi2)
{
   register double sum = 0.0;

   for(int i=0; i<np; i++)
   {
      sum += w[i]*phi1[i]*phi2[eln][i];
   }
   sum *= Jac[eln];
   return sum;
}
//-----------------------------------------------
double integr_MatMat(int np, int eln, double *w, double *Jac, double **phi1, double **phi2)
{
   register double sum = 0.0;

   for(int i=0; i<np; i++)
   {
      sum += w[i]*phi1[eln][i]*phi2[eln][i];
   }
   sum *= Jac[eln];
   return sum;
}
//===============================================
/* Assemble the M, C, L matrices and RHS 
 * f term -- For a generic element only 
 * i.e. no BCs applied
 */
//-----------------------------------------------
void assem(int Mdim, int np, int P, int eln, double *z,
 double *w, double *phi1, double *phi2, double *dphi1,
  double *dphi2, double *Jac, double **D,
  double ***M, double ***C, double ***L)
{
   /* phi1 --> phi_q
    * phi2 --> phi_p
    *    i --> q (row)
    *    j --> p (column) */
   for(int i=0; i<P+1; i++)
   {
      basis(np, P, i, z, phi1);
      diff(np, eln, phi1, D, Jac, dphi1);

      for(int j=0; j<P+1; j++)
      {
         basis(np, P, j, z, phi2);
         diff(np, eln, phi2, D, Jac, dphi2);
         M[i][j][eln] = integr(np, eln, w, Jac, phi1, phi2);
         C[i][j][eln] = integr(np, eln, w, Jac, phi1, dphi2);
         L[i][j][eln] = -integr(np, eln, w, Jac, dphi1, dphi2); // MINUS APPLIED HERE
      }
   }
}
//===============================================
/* Assemble the global grid */
//-----------------------------------------------
void xGlob(int np, int Nel, int eln, double **x, double *xG)
{
   /* Shared points will be overwritten with same value */
   for(int i=0; i<np; i++)
   {
      xG[i + eln*(np-1)] = x[eln][i];
   }
}
//===============================================
/* Mapping function: In the frequency space */
//-----------------------------------------------
void mapping_frequency(int Nel, int P, int **map)
{
   for(int i=0; i<Nel; i++)
   {
      for(int j=0; j<P+1; j++)
      {
         map[i][j] = j + i*(P+1);
      }
   }
}
//===============================================
/* Mapping function: In the physical space */
//-----------------------------------------------
void mapping_physical(int Nel, int np, int **map)
{
   for(int i=0; i<Nel; i++)
   {
      for(int j=0; j<np; j++)
      {
         map[i][j] = j + (i*(np-1));
      }
   }
}
//===============================================
/* Contribution of source term -- weak formulation */
//-----------------------------------------------
void InnerProduct_f_phi(int np, int P, int eln,
 double *phi1, double *Jac, double **f,
  double **inner_prod)
{
   for(int i=0; i<P+1; i++)
   {
      basis(np, P, i, z, phi1);

      inner_prod[eln][i] = integr_VecMat(np, eln, w, Jac, phi1, f);
   }
}
//===============================================
/* Go from weights to a given physical space */
//-----------------------------------------------
void FrequencyToPhysical(int np, int P, int eln,
   double *z, double *phi, double **u_hat, double **u_delta)
{
   /* INPUTS:
    * - np: number of points in physical space destination
    * - u_hat: weights in frequency space
    * - z: zeros (i.e. local physical space nodes)
    * OUTPUTS:
    * - u_delta: function in the desired physical space
   */

   /* Initialize as zeros */
   for(int j=0; j<np; j++)
   {
      u_delta[eln][j] = 0.0;  
   }

   /* Loop through all bases (i.e. p) */
   for(int i=0; i<(P+1); i++)
   {
      basis(np, P, i, z, phi);
      /* Loop through all points */
      for(int j=0; j<np; j++)
      {
         u_delta[eln][j] += u_hat[eln][i]*phi[j];  
      }
   }
}
//===============================================
/* Go from physical space to frequency space */
//-----------------------------------------------
void PhysicalToFrequency(int P, int eln, double *Jac, double **u_delta, double **u_hat)
{
   /* Inner product of u_delta and phi */
   InnerProduct_f_phi(np, P, eln, phi1, Jac, u_delta, DUMMY_INNER_PROD);
   /* Get the weights by solving the Ax=b problem */
   GetWeightsFromLinearAlgebra(P, eln, DUMMY_INNER_PROD, M, u_hat);
}
//===============================================
// Setup and solve Ax=b problem for weights
//-----------------------------------------------
void GetWeightsFromLinearAlgebra(int P, int eln, 
   double **ELEMENTWISE_VEC, double ***ELEMENTWISE_MAT, double **el_weights)
{
   /* INPUTS:
    * - el_weights: elementwise wise vector for storing them
    * OUTPUTS:
    * - u_delta: function in the desired physical space
   */
   /* Set up the linear algebra problem for LAPACK */
   for(int i=0; i<(P+1); i++)
   {
      DUMMY_VEC[i] = ELEMENTWISE_VEC[eln][i];
      for(int j=0; j<(P+1); j++)
      {
         DUMMY_MAT[i][j] = ELEMENTWISE_MAT[i][j][eln];
      }
   }
   /* Flatten matrix for LAPACK */
   flattenF(Mdim, Mdim, DUMMY_MAT, DUMMY_MAT_FLAT);
   /* Solve for weights using LAPACK --> Stored in 6th argument */
   dgesv_(&Mdim, &NRHS, DUMMY_MAT_FLAT, &Mdim, ipiv, DUMMY_VEC, &Mdim, &INFO);
   /* Store the solution (weights) for this element */
   for(int i=0; i<(P+1); i++)
   {
      el_weights[eln][i] = DUMMY_VEC[i];
   }
}
//===============================================