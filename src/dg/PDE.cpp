#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// polynomial library from NEKTAR
#include "polylib.h"
// My standard libraries
#include "allocatePointersLib.h"
#include "flattenLib.h"
// my libraries
#include "var.h"
#include "semlib.h"
#include "PDE.h"
#include "TimeAdvVars.h"

using namespace polylib;
using namespace std;

//===============================================
//				LAPACK FUNCTIONS
//===============================================
extern "C" {extern void dgesv_(int *, int *, double (*), int *, int [], double (*), int *, int*);}
//===============================================
void displayMainParameters()
{
	// Print main parameters
	cout << endl;
	cout << "GENERAL:" << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "~~~ Element ~~~" << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	cout << "Polynomial order: P = " << P << endl;
	cout << "Number of elements: Nel = " << Nel << endl;
	cout << "Number of points per element: np = " << np << endl;
	// cout << endl;
	// cout << "~~~ Global ~~~" << endl;
	// cout << "- - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	// cout << "Number of elements: Nel = " << Nel << endl;
	// cout << "Number of distinct points: npG = " << npG << endl;
	// cout << "Degrees of freedom (convergence): nDOF = " << nDOF << endl;
	cout << endl;
	cout << "~~~ Physics ~~~" << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	cout << "Initial time: T0 = " << T0 << endl;
	cout << "Final time: Tf = " << Tf << endl;
	cout << "CFL: " << CFL << endl;
	cout << "PDE Type: " << problem << endl;
	cout << "Subproblem: " << subproblem << endl;
	cout << "-----------------------------------------------" << endl;
	cout << endl;
	cout << "~~~ Domain ~~~" << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	cout << "Lower bound: " << xDomainBounds[0] << endl;
	cout << "Upper bound: " << xDomainBounds[1] << endl;
	cout << "-----------------------------------------------" << endl;
}
//===============================================
void updateNumPoints_ProblemSizeParameters(int P, int Nel)
{
	// Number of points + Matrix dimension
	if(P == 6)
	{
		np = P+1; // for Burgers
	}
	else if(P==7)
	{
		np = int(ceil(3.0*P/2.0 + 2.0)); // minimum number of points for incompressible NS
		// np = 17;
		// np = P+1;
	}
	else
	{
		// np = int(ceil(3.0*P/2.0 + 2.0)); // minimum number of points for incompressible NS // linear advection
		np = P+1; // For burger's
	}
	// np = int(ceil(3.0*P/2.0 + 2.0)); // minimum number of points for incompressible NS // linear advection
	/* Start fixed (cannot change - impossible) */
	npG = Nel*np-(Nel-1); // global number of points
	Mdim = P+1; // a single element
	/* End fixed */
	nDOF = Mdim*Nel; // Degrees of freedom for convergence studies
	/* Implicit */
	GlobalDim = Mdim*Nel;
}
//===============================================
double InitialCondition_Scalar(double x)
{
	register double val;
	if(problem=="LinearAdvection")
	{
		val = sin(x);
	}
	else if(problem=="NonhomogenousInviscidBurgers")
	{
		val = sin(PI*x)+0.01;
	}
	else if(problem=="HomogenousInviscidBurgers")
	{
		val = atan(100.0*(x-0.25))+atan(-200.0*(x-0.5))+atan(100.0*(x-0.7));
	}
	else if(problem=="SteadyStateLinearAdvection_Manufactured")
	{
		// val = sin(2.0*PI*x); // sine wave test
		val = 0.0;
		for(int j=0; j<arctan_NumShocks; j++)
		{
			val += atan(arctan_ShockMag[j]*(x-arctan_ShockLoc[j]));
		}
		val *= arctan_magnitude;
		val += arctan_const;
	}
	return val;
}
//===============================================
void InitialCondition(int np, int eln, double **x, double **u0)
{
	/* Elementwise */
	for(int i=0; i<np; i++)
	{
		u0[eln][i] = InitialCondition_Scalar(x[eln][i]);
	}
}
//===============================================
double SourceTerm_Scalar(double x, double t)
{
	register double val;
	if(problem=="LinearAdvection")
	{
		val = 0.0;
	}
	else if(problem=="NonhomogenousInviscidBurgers")
	{
		val = PI*cos(PI*(x-t))*(-0.99+sin(PI*(x-t)));
	}
	else if(problem=="HomogenousInviscidBurgers")
	{
		val = 0.0;
	}
	else if(problem=="SteadyStateLinearAdvection_Manufactured")
	{
		val = 0.0;
		for(int j=0; j<arctan_NumShocks; j++)
		{
			val += arctan_ShockMag[j]/(arctan_ShockMag[j]*arctan_ShockMag[j]*(x-arctan_ShockLoc[j])*(x-arctan_ShockLoc[j]) + 1.0);
		}
		val *= arctan_magnitude;
		// val = 2.0*PI*cos(2.0*PI*x); // sine wave test
	}
	return val;
}
//===============================================
double ExactSolution_Scalar(double x, double t)
{
	register double val;
	if(problem=="LinearAdvection")
	{
		val = sin(x-LinAdvSpeed*t);
	}
	else if(problem=="NonhomogenousInviscidBurgers")
	{
		val = sin(PI*(x-t))+0.01;
	}
	else if(problem=="HomogenousInviscidBurgers")
	{
		val = 0.0;
		// val = atan(100.0*(x-0.25))+atan(-200.0*(x-0.5))+atan(100.0*(x-0.7));
	}
	else if(problem=="SteadyStateLinearAdvection_Manufactured")
	{
		val = 0.0;
		for(int j=0; j<arctan_NumShocks; j++)
		{
			val += atan(arctan_ShockMag[j]*(x-arctan_ShockLoc[j]));
		}
		val *= arctan_magnitude;
		val += arctan_const;
		// val = sin(2.0*PI*x);
	}
	return val;
}
//===============================================
void SourceTerm(int np, int eln, double t, double **x, double **q_source)
{
	for(int i=0; i<np; i++)
	{
		q_source[eln][i] = SourceTerm_Scalar(x[eln][i],t);
	}
}
//===============================================
void ExactSolution(int np, int eln, double t, double **x, double **u_exact)
{
	for(int i=0; i<np; i++)
	{
		u_exact[eln][i] = ExactSolution_Scalar(x[eln][i],t);
	}
}
//===============================================
/* Flux expression */
//-----------------------------------------------
void Flux_Physical(int np, int eln, double **u, double **f)
{
	/* Elementwise */
	for(int i=0; i<np; i++)
	{
		if((problem=="LinearAdvection") || (problem=="SteadyStateLinearAdvection_Manufactured"))
		{
			f[eln][i] = LinAdvSpeed*u[eln][i];
		}
		else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
		{
			f[eln][i] = 0.5*u[eln][i]*u[eln][i];
		}
	}
}
//===============================================
/* Flux expression */
//-----------------------------------------------
double RiemannSolver(double f_L, double f_R, double u_L, double u_R)
{
	register double fstar = 0.0; // initialize
	register double maxAbsSpeed = 0.0;

	/* Lax-Friedrichs type flux */
	if((problem=="LinearAdvection") || (problem=="SteadyStateLinearAdvection_Manufactured"))
	{
		maxAbsSpeed = abs(LinAdvSpeed);
	}
	else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
	{
		if(abs(u_R)>abs(u_L))
		{
			maxAbsSpeed = abs(u_R);
		}
		else
		{
			maxAbsSpeed = abs(u_L);
		}
	}
	fstar = 0.5*(f_L+f_R) - 0.5*maxAbsSpeed*(u_R-u_L);

	/* Upwind type flux */
	// if(problem=="LinearAdvection")
	// {
	// 	fstar = f_L;
	// }
	// else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
	// {
		/* Upwind 0 */
		// maxAbsSpeed = (f_R - f_L)/(u_R - u_L);
		// fstar = 0.5*(f_L+f_R) - 0.5*abs(maxAbsSpeed)*(u_R-u_L);
		/* Upwind 1 */
		// if(maxAbsSpeed < 0.0)
		// {
		// 	fstar = f_R;
		// }
		// else
		// {
		// 	fstar = f_L;
		// }
		/* Upwind 2 */
		// if(u_L < 0.0)
		// {
		// 	fstar = f_R;
		// }
		// else
		// {
		// 	fstar = f_L;
		// }
		/* Upwind 3 */
		// maxAbsSpeed = 0.5*(u_R+u_L); // test
		// if(maxAbsSpeed < 0.0)
		// {
		// 	fstar = f_R;
		// }
		// else
		// {
		// 	fstar = f_L;
		// }
		/* ROE type flux */
		// fstar = 0.5*(f_L+f_R) - 0.5*abs(u_R + u_L)*(u_R-u_L);
		/* ECON flux */
		// fstar = 0.5*(f_L+f_R) - (1.0/12.0)*(u_R-u_L)*(u_R-u_L);
	// }
	
	return fstar;
}
//===============================================
/* Flux expression */
//-----------------------------------------------
void Flux_Numerical(int np, int eln, int Nel, double TIME, double **u, double **f, double **f_star)
{
	register double f_L, f_R, u_L, u_R;
	
	/* Elementwise */

	// --------------------------------
	/* At the LEFT FACE of element */
	// --------------------------------
	if(eln==0)
	{
		/* Left boundary element */
		// Dirichlet BC at left boundary (x=xL)
		// - Solution
		if(problem=="HomogenousInviscidBurgers")
		{
			// u_L = u[Nel-1][np-1]; // periodic BCs
			u_L = ExactSolution_Scalar(xDomainBounds[0],TIME); // Dirichlet test
		}
		else
		{
			u_L = ExactSolution_Scalar(xDomainBounds[0],TIME);
		}
		// - Flux
		if((problem=="LinearAdvection") || (problem=="SteadyStateLinearAdvection_Manufactured"))
		{
			f_L = LinAdvSpeed*u_L; // ghost element
		}
		else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
		{
			f_L = 0.5*u_L*u_L;
		}
	}
	else
	{
		/* Interior elements */
		u_L = u[eln-1][np-1];	
		f_L = f[eln-1][np-1];
	}
	u_R = u[eln][0];
	f_R = f[eln][0];
	// f_star[eln][0] = 0.5*(f_L+f_R) + 0.5*abs(LinAdvSpeed)*(1.0-alpha)*(-1.0*f_L+1.0*f_R);
	f_star[eln][0] = RiemannSolver(f_L,f_R,u_L,u_R);
	
	// --------------------------------
	/* At the RIGHT FACE of element */
	// --------------------------------
	if(eln==(Nel-1))
	{
		/* Right boundary element */
		// Dirichlet BC at right boundary (x=xR)
		// - Solution
		if(problem=="HomogenousInviscidBurgers")
		{
			// u_R = u[0][0]; // periodic BCs
			u_R = ExactSolution_Scalar(xDomainBounds[1],TIME); // Dirichlet, test
		}
		else
		{
			u_R = ExactSolution_Scalar(xDomainBounds[1],TIME);
		}
		// - Flux
		if((problem=="LinearAdvection") || (problem=="SteadyStateLinearAdvection_Manufactured"))
		{
			f_R = LinAdvSpeed*u_R; // ghost element
		}
		else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
		{
			f_R = 0.5*u_R*u_R; // ghost element
		}
	}
	else
	{
		/* Interior elements */
		// At the right boundary:
		u_R = u[eln+1][0];
		f_R = f[eln+1][0];
	}
	u_L = u[eln][np-1];
	f_L = f[eln][np-1];
	// f_star[eln][1] = 0.5*(f_L+f_R) + 0.5*abs(LinAdvSpeed)*(1.0-alpha)*(-1.0*f_L+1.0*f_R);
	f_star[eln][1] = RiemannSolver(f_L,f_R,u_L,u_R);
}
//===============================================
/* Compute RHS of strong form PDE as a fxn of u */
//-----------------------------------------------
void UpdateWeights(double TIME, double TIME_STEP, double **u_hat, double **u_hat_new)
{
	/* Update physical fluxes for all elements first
	 * because numerical flux depends on
	 * neighboring elements
	*/
	for(int eln=0; eln<Nel; eln++)
	{
		/* Transform to the physical space */	
		FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_dummy);
		/* Compute physical flux (f) */
		Flux_Physical(np, eln, u_dummy, f);
		/* Compute source term */
		SourceTerm(np, eln, TIME, x, q_source);
		/* Compute inner product of source term and phi */
		InnerProduct_f_phi(np, P, eln, phi1, Jac, q_source, q_phi_innerproduct);
	}

	/* Obtain the new weights */
	for(int eln=0; eln<Nel; eln++)
	{
		/* Compute numerical flux (f_star) */
		Flux_Numerical(np, eln, Nel, TIME, u_dummy, f, f_star);
		/* Transform physical flux to freq. space */
		PhysicalToFrequency(P, eln, Jac, f, f_hat);
		/* Compute RHS of the PDE: 
		 * - NOTE: Strong form is used here
		*/
		for(int i=0; i<(P+1); i++)
		{
			RHS_VEC[eln][i] = q_phi_innerproduct[eln][i]
							+(f[eln][np-1]-f_star[eln][1])*phi_R[eln][i]
							-(f[eln][0]-f_star[eln][0])*phi_L[eln][i];
			
			if(problem!="SteadyStateLinearAdvection_Manufactured")
			{
				for(int j=0; j<(P+1); j++)
				{
					RHS_VEC[eln][i] += -C[i][j][eln]*f_hat[eln][j];
				}
				RHS_VEC[eln][i] *= TIME_STEP;
			}
		}
		/* Get weights by solving Ax=b problem */
		if(problem!="SteadyStateLinearAdvection_Manufactured")
		{
			GetWeightsFromLinearAlgebra(P, eln, RHS_VEC, M, u_hat_new);
		}
		else
		{ // remove this and all the rest of the stuff for SS LA method 1
			for(int i=0; i<(P+1); i++)
			{
				for(int j=0; j<(P+1); j++)
				{
					L[i][j][eln] = C[i][j][eln] + M[i][j][eln];
				}
			}
			/* Case: SteadyStateLinearAdvection_Manufactured */
			GetWeightsFromLinearAlgebra(P, eln, RHS_VEC, M, u_hat_new);
		}
	}
}
// //===============================================
// /* TIME ADVANCEMENT: Explicit */
// //-----------------------------------------------
void ExplicitTimeAdv(double TIME, double TIME_STEP, double **u_hat)
{
	UpdateWeights(TIME, TIME_STEP, u_hat, K1);
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			u_hat[eln][i] = u_hat[eln][i] + K1[eln][i];
		}
	}
}
//===============================================
/* TIME ADVANCEMENT: RK4 */
//-----------------------------------------------
void RK4(double TIME, double TIME_STEP, double **u_hat)
{
	// STEP 1:
	UpdateWeights(TIME, TIME_STEP, u_hat, K1);
	// STEP 2:
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			RK_INPUT_VEC[eln][i] = u_hat[eln][i] + 0.5*K1[eln][i];
		}
	}
	UpdateWeights((TIME+0.5*TIME_STEP), TIME_STEP, RK_INPUT_VEC, K2);
	// STEP 3:
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			RK_INPUT_VEC[eln][i] = u_hat[eln][i] + 0.5*K2[eln][i];
		}
	}
	UpdateWeights((TIME+0.5*TIME_STEP), TIME_STEP, RK_INPUT_VEC, K3);
	// STEP 4:
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			RK_INPUT_VEC[eln][i] = u_hat[eln][i] + K3[eln][i];
		}
	}
	UpdateWeights((TIME+TIME_STEP), TIME_STEP, RK_INPUT_VEC, K4);
	// FINAL STEP: -- DELTA_T is applied inside UpdateWeights()
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			u_hat[eln][i] = u_hat[eln][i] + (1.0/6.0)*(K1[eln][i]+2.0*K2[eln][i]+2.0*K3[eln][i]+K4[eln][i]);
		}
	}
}
//===============================================
/* Update time step */
//-----------------------------------------------
double getTimeStep_Burgers(double **u_hat)
{
	register double dt, abs_u_max = 0.0;

	for(int eln=0; eln<Nel; eln++)
	{
		/* Transform to the physical space */	
		FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_dummy); // temporary

		for(int i=0; i<np; i++)
		{
			if(abs(u_dummy[eln][i]) > abs_u_max)
			{
				abs_u_max = abs(u_dummy[eln][i]);
			}
		}
	}
	dt = CFL*DELTA_X_MIN/abs_u_max;
	return dt;
}
//===============================================
/* Implicit stuff below */
//-----------------------------------------------
void ElementwiseToGlobal_Matrix(int P, int Nel, int **map_frequency, double ***ElementwiseMatrix, double **GlobalMatrix)
{
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			for(int j=0; j<(P+1); j++)
			{
				GlobalMatrix[map_frequency[eln][i]][map_frequency[eln][j]] = ElementwiseMatrix[i][j][eln];
			}
		}
	}
}
//===============================================
void ElementwiseToGlobal_Vector(int P, int Nel, int **map_frequency, double **ElementwiseVector, double *GlobalVector)
{
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			GlobalVector[map_frequency[eln][i]] = ElementwiseVector[eln][i];
		}
	}
}
//===============================================
void GlobalToElementwise_Vector(int P, int Nel, int **map_frequency, double *GlobalVector, double **ElementwiseVector)
{
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			ElementwiseVector[eln][i] = GlobalVector[map_frequency[eln][i]];
		}
	}
}
//===============================================
void assemGlobalVectorRHS(double *GlobalVector)
{
	/* Case: SteadyStateLinearAdvection_Manufactured */
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			GlobalVector[map_frequency[eln][i]] = q_phi_innerproduct_global[map_frequency[eln][i]];
			if(eln==0)
			{
				GlobalVector[map_frequency[eln][i]] += ExactSolution_Scalar(xDomainBounds[0],0.0)*phi_L[eln][i];
			}
		}
	}	
}
//===============================================
void assemGlobalMatrixLHS(double **GlobalMatrix)
{
	/* Case: SteadyStateLinearAdvection_Manufactured */
	for(int i=0; i<GlobalDim; i++)
	{
		for(int j=0; j<GlobalDim; j++)
		{
			GlobalMatrix[i][j] = C_global[i][j] + PHI_MAT_global[i][j];
		}
	}
}
//===============================================
void compute_PHI_MATS()
{
	// Compute the PHI_LL_MAT and PHI_LR_MAT
	/* Case: SteadyStateLinearAdvection_Manufactured */
	for(int i=0; i<(P+1); i++)
	{
		for(int j=0; j<(P+1); j++)
		{
			PHI_LL_MAT[i][j] = phi_L[0][i]*phi_L[0][j];
			PHI_LR_MAT[i][j] = phi_L[0][i]*phi_R[0][j];
		}	
	}	
}
//===============================================
void assemGlobalPhiMatrix()
{
	// Compute PHI_MAT_global Matrix
	/* Case: SteadyStateLinearAdvection_Manufactured */
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<(P+1); i++)
		{
			for(int j=0; j<(P+1); j++)
			{
				PHI_MAT_global[map_frequency[eln][i]][map_frequency[eln][j]] = PHI_LL_MAT[i][j];
				if(eln > 0)
				{
					PHI_MAT_global[map_frequency[eln][i]][map_frequency[eln-1][j]] = -PHI_LR_MAT[i][j];
				}
			}
		}
	}
}
//===============================================

//===============================================
/* Setup the Time Adv Parameters */
//-----------------------------------------------
// void adjustTimeStep(double CFL_input, double DELTA_X_MIN)
// {
// 	// DECLARE IN HEADER
// 	TIME_STEP = (CFL_input*DELTA_X_MIN)/LinAdvSpeed;

// 	register double t = T0; // dummy var
// 	NOUT = 0; // initialize

// 	while(t<Tf)
// 	{
// 		t += TIME_STEP;
// 		NOUT += 1;
// 	}
// 	// cout << "CFL: " << CFL << endl;
// 	// cout << "NOUT: " << NOUT << endl;
// 	// cout << "t: " << t << endl;
// 	// cout << "TIME_STEP: " << TIME_STEP << endl;

// 	/* Correct time step for exactly reaching Tf */
// 	// - This will be slightly smaller
// 	TIME_STEP = (Tf-T0)/double(NOUT);
// 	CFL = 0.5*TIME_STEP/DELTA_X_MIN; // update CFL

// 	t = T0;
// 	NOUT = 0;
// 	while(t<Tf)
// 	{
// 		t += TIME_STEP;
// 		NOUT += 1;
// 	}
// 	// cout << "CFL: " << CFL << endl;
// 	// cout << "NOUT: " << NOUT << endl;
// 	// cout << "t: " << t << endl;
// 	// cout << "TIME_STEP: " << TIME_STEP << endl;
// }
// //-----------------------------------------------