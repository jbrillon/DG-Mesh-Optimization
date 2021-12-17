SetupFile="./SolverSetupFiles/MDO.setup"
ParametersFile="./parameters/arctangent.parameters"
P="5"
# get .exe after building code
cp ../../1D_DG_solver_01.exe ./
# run
rm h_type_extension.txt
touch h_type_extension.txt
# - h convergence
# Type of convergence:
# - (0==p, 1==h)
echo "Bash version ${BASH_VERSION}..."
# inputs (ordered): ConvStudyFlag, P, Nel, ConvStudyType (0==p, 1==h)
for i in $(seq 2 4 100)
do
	# echo -e "1\n4\n$i\n1\n" | ./1D_DG_solver_01.exe
	echo -e "1\n$P\n$i\n1\n$SetupFile\n$ParametersFile\n" | ./1D_DG_solver_01.exe
done
for i in $(seq 100 10 200)
do
	# echo -e "1\n4\n$i\n1\n" | ./1D_DG_solver_01.exe
	echo -e "1\n$P\n$i\n1\n$SetupFile\n$ParametersFile\n" | ./1D_DG_solver_01.exe
done
for i in $(seq 220 20 300)
do
	# echo -e "1\n4\n$i\n1\n" | ./1D_DG_solver_01.exe
	echo -e "1\n$P\n$i\n1\n$SetupFile\n$ParametersFile\n" | ./1D_DG_solver_01.exe
done
for i in $(seq 350 50 500)
do
	# echo -e "1\n4\n$i\n1\n" | ./1D_DG_solver_01.exe
	echo -e "1\n$P\n$i\n1\n$SetupFile\n$ParametersFile\n" | ./1D_DG_solver_01.exe
done