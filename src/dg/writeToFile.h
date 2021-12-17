#ifndef WRITETOFILE_H
#define WRITETOFILE_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>

void writeIntToFile(int, std::string, std::string, std::string);
void writeDoubleToFile(double, std::string, std::string, std::string);
void writeStringArrayToFile(int, std::string*, std::string, std::string, std::string);
void write2ArraysToFile(int, double*, double*, std::string, std::string, std::string);
void writeArrayToFile(int, double*, std::string, std::string, std::string);
void writeMatrixToFile(int, int, double**, std::string, std::string, std::string);
void writeTransientSolutionToFile(int, int, double **, 
	std::string, std::string, std::string, int, double);

#endif