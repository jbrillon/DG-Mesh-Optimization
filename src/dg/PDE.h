#ifndef PDE_H
#define PDE_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dglib.h"

void displayMainParameters();
void updateNumPoints_ProblemSizeParameters(int, int);
double InitialCondition_Scalar(double);
void InitialCondition(int, int, double **, double **);
double SourceTerm_Scalar(double, double);
double ExactSolution_Scalar(double, double);
void SourceTerm(int, int, double, double **, double **);
void ExactSolution(int, int, double, double **, double **);
void Flux_Physical(int, int, double **, double **);
double RiemannSolver(double, double, double, double);
void Flux_Numerical(int, int, int, double, double **, double **, double **);
void UpdateWeights(double, double, double**, double**);
void ExplicitTimeAdv(double, double, double **);
void RK4(double, double, double **);
double getTimeStep_Burgers(double **);
void ElementwiseToGlobal_Matrix(int, int, int **, double ***, double **);
void ElementwiseToGlobal_Vector(int, int, int **, double **, double *);
void GlobalToElementwise_Vector(int, int, int **, double *, double **);
void assemGlobalVectorRHS(double *);
void assemGlobalMatrixLHS(double **);
void compute_PHI_MATS();
void assemGlobalPhiMatrix();
// void adjustTimeStep(double, double);

#endif