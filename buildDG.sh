cd ./src/dg/
# includes
g++ -c polylib.c
g++ -c allocatePointersLib.cpp
g++ -c writeToFile.cpp
g++ -c flattenLib.cpp
g++ -c strManipLib.cpp
g++ -c var.cpp
g++ -c allocateVar.cpp
g++ -c dglib.cpp
g++ -c ErrorEstimates.cpp
g++ -c PDE.cpp
g++ -c readFilesDG.cpp
g++ -c initZeros.cpp
# main
g++ -c main.cpp
g++ -o DG_solver.exe main.o initZeros.o readFilesDG.o PDE.o ErrorEstimates.o dglib.o allocateVar.o var.o strManipLib.o flattenLib.o writeToFile.o allocatePointersLib.o polylib.o -llapack -lblas -framework Accelerate
# remove all .o files
rm *.o
# move .exe to root
mv DG_solver.exe ../../