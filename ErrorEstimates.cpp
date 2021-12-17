#include <iostream>
#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <math.h>
#include <iomanip>
// polynomial library from NEKTAR
#include "polylib.h"
// My standard libraires:
#include "allocatePointersLib.h"
#include "writeToFile.h"
#include "flattenLib.h"
// My SEM libraries
#include "var.h"
#include "semlib.h"
#include "ErrorEstimates.h"

using namespace std;

//-----------------------------------------------
double L2_NORM(double **u_delta,double **u_exact)
{
	register double sum = 0.0;
	/* Compute the difference between solutions */
	for(int eln=0; eln<Nel; eln++)
	{
		for(int i=0; i<np; i++)
		{
			u_dummy[eln][i]	= u_delta[eln][i] - u_exact[eln][i];
		}	
	}
	/* Sum the integral of difference squared */
	for(int eln=0; eln<Nel; eln++)
	{
		sum += integr_MatMat(np, eln, w, Jac, u_dummy, u_dummy);
	}
	/* Take square root to obtain L2-norm */
	sum = pow(sum,0.5);
	return sum;
}
//-----------------------------------------------