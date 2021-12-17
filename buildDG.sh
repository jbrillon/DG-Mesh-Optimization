# includes
g++ -c polylib.c
g++ -c allocatePointersLib.cpp
g++ -c writeToFile.cpp
g++ -c flattenLib.cpp
g++ -c strManipLib.cpp
g++ -c var.cpp
g++ -c allocateVar.cpp
g++ -c semlib.cpp
g++ -c ErrorEstimates.cpp
g++ -c PDE.cpp
g++ -c readFilesDG.cpp
g++ -c initZeros.cpp
# main
g++ -c 1D_DG_solver_01.cpp
g++ -o 1D_DG_solver_01.exe 1D_DG_solver_01.o initZeros.o readFilesDG.o PDE.o ErrorEstimates.o semlib.o allocateVar.o var.o strManipLib.o flattenLib.o writeToFile.o allocatePointersLib.o polylib.o -llapack -lblas -framework Accelerate