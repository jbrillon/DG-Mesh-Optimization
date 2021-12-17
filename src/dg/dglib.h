#ifndef DGLIB_H
#define DGLIB_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void elementBounds(int, double*, double*);
void chi(int, int, double*, double*, double**, double*);
void basis(int, int, int, double*, double*);
void basis_vec_xLR(int, int, double*, double*, double**, double**);
void diff(int, int, double*, double**, double*, double*);
double integr(int, int, double*, double*, double*, double*);
double integr_VecMat(int, int, double*, double*, double*, double**);
double integr_MatMat(int, int, double*, double*, double**, double**);
void assem(int, int, int, int, 
	double *, double *, double *, 
	double *, double *, double *, double *, double **,
	double ***, double ***, double ***);
void xGlob(int, int, int, double**, double*);
void mapping_frequency(int, int, int**);
void mapping_physical(int, int, int**);
void InnerProduct_f_phi(int, int, int, double *, double *, double **, double **);
void FrequencyToPhysical(int, int, int, double*, double*, double**, double**);
void PhysicalToFrequency(int, int, double*, double**, double**);
void GetWeightsFromLinearAlgebra(int, int, double**, double***, double **);

#endif