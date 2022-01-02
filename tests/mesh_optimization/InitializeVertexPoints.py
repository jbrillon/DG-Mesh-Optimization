# Import libraries
import numpy as np # NumPy: contains basic numerical routines

SetupFile="./mesh_optimization.setup"

# def readSetupFile():
file = open(SetupFile,"r")
line = file.readline().replace('\n','').strip()
while(line != 'START'):
	line = file.readline().replace('\n','').strip()
line = file.readline().replace('\n','').strip()
while(line != 'END'):
	line = file.readline().replace('\n','').strip()
	if(line == "NUMBER OF ELEMENTS"):
		line = file.readline().replace('\n','').strip()
		Nel = int(line)
	elif(line == "LEFT BOUNDARY X-VALUE"):
		line = file.readline().replace('\n','').strip()
		xL = np.float64(line)
	elif(line == "RIGHT BOUNDARY X-VALUE"):
		line = file.readline().replace('\n','').strip()
		xR = np.float64(line)

print("Nel")
print(Nel)
print("xL")
print(xL)
print("xR")
print(xR)

xgrid = np.linspace(xL,xR,Nel+1)
xvert = xgrid[1:Nel]

def writeVerticesToFile(vertices,filename,data_fileType):
	file = open(filename+"."+data_fileType,"w")
	file.write("START")
	file.write("\n")
	for v in vertices:
		wstr = "%.16e" % v
		file.write(wstr)
		file.write("\n")
	file.write("END")

writeVerticesToFile(xvert,"vertices","txt")

import subprocess # as sp
print('=====================================================')

SetupFile = "./mesh_optimization.setup"
ParametersFile = "../../parameters/arctangent.parameters"
proc = subprocess.Popen("./1D_DG_solver_01.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
proc.communicate("0\n0\n"+SetupFile+"\n"+"./vertices.txt\n"+ParametersFile+"\n")
# proc.communicate("0\n0\n"+SetupFile+"\n"+ParametersFile+"\n")
proc.wait()
L2error = np.loadtxt("./Data/L2error.txt")
print("From Python:")
print(L2error)

err_file = open("L2error_python.txt",'w')

wstr = "%.16e\n" % L2error
err_file.write(wstr)

err_file.close()

print('=====================================================')