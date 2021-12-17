#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <math.h>
#include <iomanip> 
#include <stdio.h>
// My libraires:
#include "var.h"
#include "readFilesDG.h"
#include "strManipLib.h"
#include "allocatePointersLib.h"

/* Namespaces */
using namespace std;

//===============================================
void readSetupFile(string setupFilename)
{
    string line;

	ifstream SETUPF (setupFilename);
	
	getline(SETUPF, line);
	while(line != "START")
	{
		getline(SETUPF, line);	
	}
	getline(SETUPF, line);

	while(line != "END")
	{
		/* Read P and Nel only if 
		   a convergence study is NOT being done */
		if((line == "POLYNOMIAL ORDER") && (ConvStudyFlag == 0))
		{
			getline(SETUPF, line);
			P = (int)stof(line);
		}
		else if((line == "NUMBER OF ELEMENTS") && (ConvStudyFlag == 0))
		{
			getline(SETUPF, line);
			Nel = (int)stof(line);
		}
		else if(line == "CFL")
		{
			getline(SETUPF, line);
			CFL = stod(line);
		}
		else if(line == "PDE TYPE")
		{
			getline(SETUPF, line);
			problem = line;
		}
		else if(line == "SUB-PROBLEM")
		{
			getline(SETUPF, line);
			subproblem = line;
		}
		else if(line == "INITIAL TIME")
		{
			getline(SETUPF, line);
			T0 = stod(line);
		}
		else if(line == "FINAL TIME")
		{
			getline(SETUPF, line);
			Tf = stod(line);
		}
		else if(line == "LEFT BOUNDARY X-VALUE")
		{
			getline(SETUPF, line);
			xDomainBounds[0] = stod(line);
		}
		else if(line == "RIGHT BOUNDARY X-VALUE")
		{
			getline(SETUPF, line);
			xDomainBounds[1] = stod(line);
		}
		getline(SETUPF, line);
	}
	SETUPF.close();
}
//===============================================
void readElementVertices(string elvertFilename)
{
    string line;

	ifstream ELVERTF (elvertFilename);
	
	getline(ELVERTF, line);
	while(line != "START")
	{
		getline(ELVERTF, line);	
	}
	getline(ELVERTF, line);
	elBounds[0] = xDomainBounds[0];
	int i = 1;
	while(line != "END")
	{
		if(i>Nel)
		{
			cout << "ERROR: Too many vertices listed in input file." << endl;
			break;
		}
		elBounds[i] = stod(line);
		i += 1;
		getline(ELVERTF, line);
	}

	elBounds[Nel] = xDomainBounds[1];

	ELVERTF.close();
}
//===============================================
void readArctanParameters(string setupFilename)
{
    string line;
    string::size_type sz1;

	ifstream SETUPF (setupFilename);
	
	getline(SETUPF, line);
	while(line != "START")
	{
		getline(SETUPF, line);	
	}
	getline(SETUPF, line);

	while(line != "END")
	{
		/* Read P and Nel only if 
		   a convergence study is NOT being done */
		if(line == "A")
		{
			getline(SETUPF, line);
			arctan_const = stod(line);
		}
		else if(line == "MAGNITUDE")
		{
			getline(SETUPF, line);
			arctan_magnitude = stod(line);
		}
		else if(line == "NUMBER OF SHOCKS")
		{
			getline(SETUPF, line);
			arctan_NumShocks = (int)stof(line);
			arctan_ShockMag = dvector(arctan_NumShocks);
			arctan_ShockLoc = dvector(arctan_NumShocks);
		}
		else if(line == "SHOCK STEEPNESS")
		{
			getline(SETUPF, line);
			sz1 = 0;
			for(int i=0; i<arctan_NumShocks; i++)
			{
				arctan_ShockMag[i] = stod(line,&sz1);
				line = line.substr(sz1);
			}
		}
		else if(line == "SHOCK LOCATIONS")
		{
			getline(SETUPF, line);
			sz1 = 0;
			for(int i=0; i<arctan_NumShocks; i++)
			{
				arctan_ShockLoc[i] = stod(line,&sz1);
				line = line.substr(sz1);
			}
		}
		getline(SETUPF, line);
	}
	SETUPF.close();
}
//===============================================