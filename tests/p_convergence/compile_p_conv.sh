# get .exe after building code
cp ../../1D_DG_solver_01.exe ./
# run
rm p_type_extension.txt
touch p_type_extension.txt
# - p convergence
# Type of convergence:
# - (0==p, 1==h)
# ConvStudyType = 0
# Nel = 2
echo "Bash version ${BASH_VERSION}..."
# inputs (ordered): ConvStudyFlag, P, Nel, ConvStudyType (0==p, 1==h)
for i in $(seq 1 2 30)
do
	echo -e "1\n$i\n2\n0\n" | ./1D_DG_solver_01.exe
done