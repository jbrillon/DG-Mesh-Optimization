#include <iostream>
#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include "writeToFile.h"
using namespace std;
//===============================================
void writeIntToFile(int flag, string fileName,
 string fileType, string directory)
{
    /* 
    Description: Writes an integer to a file 
    
    Inputs:
        flag = integer value integer
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        myfile << flag << "\n";
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void writeDoubleToFile(double val, string fileName,
 string fileType, string directory)
{
    /* 
    Description: Writes a double to a file 
    
    Inputs:
        val = double value
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        myfile << setprecision(16) << scientific << val << "\n";
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void writeStringArrayToFile(int size, string *x,
 string fileName, string fileType, string directory)
{
    /* 
    Description: Writes an string array to a file 
    
    Inputs:
        size = size of array
        x = string array
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        for(int count = 0; count < size; count ++)
        {
            myfile << x[count] << "\n";
        }

        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void write2ArraysToFile(int size, double *x, double *y,
 string fileName, string fileType, string directory)
{
    /* 
    Description: Writes an array to a file 
    
    Inputs:
        size = size of array
        x = array
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        for(int count = 0; count < size; count ++)
        {
            myfile << setprecision(16) << x[count] << "\t" << y[count] << "\n";
        }
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void writeArrayToFile(int size, double *x,
 string fileName, string fileType, string directory)
{
    /* 
    Description: Writes an array to a file 
    
    Inputs:
        size = size of array
        x = array
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        for(int count = 0; count < size; count ++)
        {
            myfile << setprecision(16) << x[count] << "\n";
        }
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void writeMatrixToFile(int nRow, int nCol, double **A,
 string fileName, string fileType, string directory)
{
    /* 
    Description: Writes an array to a file 
    
    Inputs:
        nRow = number of rows
        nCol = number of columns
        A  = matrix
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    ofstream myfile (directory + "/" + fileName + "." + fileType);
    
    if (myfile.is_open())
    {
        for(int i=0; i<nRow; i++)
        {
            for(int j=0; j<nCol; j++)
            {
                myfile << A[i][j] << "\t";
            }
            myfile << "\n";
        }
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================
void writeTransientSolutionToFile(int nRow, int nCol, double **A,
 string fileName, string fileType, string directory,
 int iout, double t)
{
    /* 
    Description: Writes an array to a file 
    
    Inputs:
        nRow = number of rows
        nCol = number of columns
        A  = matrix
        fileName = name of file
        fileType = file extension (txt, dat, ...)
    */

    register string output_number;
    
    if(iout < 10)
    {
        output_number = "000" + to_string(iout);
    }
    else if(iout < 100)
    {
        output_number = "00" + to_string(iout);
    }
    else if(iout < 1000)
    {
        output_number = "0" + to_string(iout);
    }

    ofstream myfile (directory + "/" + fileName + "/" + output_number + "." + fileType);
    
    if (myfile.is_open())
    {
        myfile << t << "\n";
        for(int i=0; i<nRow; i++)
        {
            for(int j=0; j<nCol; j++)
            {
                myfile << A[i][j] << "\t";
            }
            myfile << "\n";
        }
        
        myfile.close();
    }
    else cout << "Unable to open file" << endl;
}
//===============================================