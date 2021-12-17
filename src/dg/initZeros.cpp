#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// My standard libraries
#include "initZeros.h"
//===============================================
using namespace std;
//===============================================
void initZeros_dmatrix(int nRow, int nCol, double **MAT)
{
    for(int i=0; i<nRow; i++)
    {
        for(int j=0; j<nCol; j++)
        {
            MAT[i][j] = 0.0;
        }
    }   
}
//===============================================
void initZeros_dvector(int nRow, double *VEC)
{
    for(int i=0; i<nRow; i++)
    {
        VEC[i] = 0.0;
    }
}
//===============================================