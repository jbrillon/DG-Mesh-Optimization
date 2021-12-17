#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flattenLib.h"
//===============================================
void flattenF(int nRow, int nCol, double **M, double *M_flat)
{
	/* Flatten matrix with FORTRAN-ordering */
	for(int j=0; j<nCol; j++)
	{
		for(int i=0; i<nRow; i++)
		{
			M_flat[j*nRow + i] = M[i][j];
		}
	}
}
//===============================================
void unflattenF(int nRow, int nCol, double **M, double *M_flat)
{
	/* Unflatten matrix with FORTRAN-ordering */
	for(int j=0; j<nCol; j++)
	{
		for(int i=0; i<nRow; i++)
		{
			M[i][j] = M_flat[j*nRow + i];
		}
	}
}
//===============================================
void flattenC(int nRow, int nCol, double **M, double *M_flat)
{
	/* Flatten matrix with C-ordering */
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
		{
			M_flat[i*nCol + j] = M[i][j];
		}
	}
}
//===============================================
void unflattenC(int nRow, int nCol, double **M, double *M_flat)
{
	/* Unflatten matrix with C-ordering */
	for(int i=0; i<nRow; i++)
	{
		for(int j=0; j<nCol; j++)
		{
			M[i][j] = M_flat[i*nCol + j];
		}
	}
}
//===============================================