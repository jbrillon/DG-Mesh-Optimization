# Import libraries
import numpy as np
from numpy import linalg as LA
import subprocess
#=====================================================
def readSetupFile(SetupFileName):
    global Nel, xL, xR
    file = open(SetupFileName,"r")
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
#=====================================================
def writeHessianToFile(H,k):
    global H_store_directory
    name = H_store_directory + "H_iter_%i" % k
    # write hessian to txt file
    np.savetxt(name+".txt",H,fmt='%.16e')
#=====================================================
def writeVerticesToFile(vertices,file):
    file = open(file,"w")
    file.write("START")
    file.write("\n")
    for v in vertices:
        wstr = "%.16e" % v
        file.write(wstr)
        file.write("\n")
    file.write("END")
    file.close()
#=====================================================
def appendDoubleToFile(val,file_var):
    wstr = "%.16e\n" % val
    file_var.write(wstr)
#=====================================================
def appendIntToFile(val,file_var):
    wstr = "%i\n" % val
    file_var.write(wstr)
#=====================================================
def appendVerticesToFile(x_vertices):
    global vert_file, data_dir
    vert_file = open(data_dir+"vertices_MDO_store.txt",'a')
    for v in x_vertices:
        wstr = "%.16e" % v
        vert_file.write(wstr)
        wstr = " "
        vert_file.write(wstr)
    wstr = "\n"
    vert_file.write(wstr)
    # Close the files for storing stuff
    vert_file.close()
#=====================================================
def appendToEachFile(err_val,L2error_val,iter_val,x_vertices):
    global L2error_file, vert_file, err_file, iter_file, data_dir
    L2error_file = open(data_dir+"L2error_MDO_store.txt",'a')
    err_file = open(data_dir+"err_MDO_store.txt",'a')
    iter_file = open(data_dir+"iter_MDO_store.txt",'a')

    appendVerticesToFile(x_vertices)
    wstr = "%i\n" % iter_val
    iter_file.write(wstr)
    # appendIntToFile(iter_val,iter_file)
    # appendDoubleToFile(err_val,err_file)
    wstr = "%.16e\n" % err_val
    err_file.write(wstr)
    # appendDoubleToFile(L2error_val,L2error_file)
    wstr = "%.16e\n" % L2error_val
    L2error_file.write(wstr)

    # Close the files for storing stuff
    L2error_file.close()
    err_file.close()
    iter_file.close()
#=====================================================
def obj_func(xvert):
    global SetupFile, ParametersFile, VerticesFile, VerticesFile_sh

    # Write the vertices to the file for C++ to read
    writeVerticesToFile(xvert,VerticesFile)
    # Run DG code with xvert as input
    proc = subprocess.Popen("./1D_DG_solver_01.exe", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.communicate("0\n0\n"+SetupFile+"\n"+VerticesFile_sh+"\n"+ParametersFile+"\n")
    # Get the output of the Python code
    L2error = np.loadtxt("./Data/L2error.txt")

    return L2error
#=====================================================
#               Numerical gradient
#-----------------------------------------------------
def gradient(xvert):
    global nVert
    h = 1.0e-5 # perturbation size
    gradient_val = np.zeros(nVert)
    # Numerical scheme: FD2
    # - Finite difference 2nd order
    for i in range(0,nVert):
        xvert_perturb_plus = 1.0*xvert
        xvert_perturb_minus = 1.0*xvert
        xvert_perturb_plus[i] += h
        xvert_perturb_minus[i] -= h
        f_perturb_plus = obj_func(xvert_perturb_plus)
        f_perturb_minus = obj_func(xvert_perturb_minus)
        gradient_val[i] = (f_perturb_plus-f_perturb_minus)/(2.0*h)
    return gradient_val
#=====================================================
#           Quasi Newton BFGS Solver
#-----------------------------------------------------
def QuasiNewtonMethod(X0,H):
    global sol_store, err_store, iter_store, numvars, maxIter, print_rate, tol
    global L2error_file, vert_file, err_file, iter_file, data_dir
    global restart_switch, iter_restart
    X = 1.0*X0
    g0 = gradient(X)
    err = LA.norm(g0)
    # initialize iteration counter
    if(restart_switch==0):
        k = 0
    else:
        k = 1*iter_restart

    # STEP SIZE -- USER INPUT
    alpha = 0.025 # restart_level == 0
    # alpha = 0.025 # restart_level == 1
    # alpha = 0.05 # too high

    L2error = obj_func(X)
    print("k: %i\t err: %.3e\t L2error: %.6e" % (k,err,L2error))

    if(restart_switch==0):
        appendToEachFile(err,L2error,k,X)

    while((err > tol) and ((k+1)<maxIter)):
        
        if(LA.norm(g0)==0.0):
            break
        else:
            d = -np.dot(H,g0)

        writeHessianToFile(H,k)

        delta_X = alpha*d
        X += delta_X

        g1 = gradient(X)
        delta_g = g1 - g0

        # DFP:
        # dH = (np.outer(delta_X,delta_X)/np.dot(delta_X,delta_g))-np.outer(np.dot(H,delta_g),np.dot(H,delta_g))/np.dot(delta_g,np.dot(H,delta_g))
        # BFGS:
        dH = (1.0+np.dot(delta_g,np.dot(H,delta_g))/np.dot(delta_g,delta_X))*(np.outer(delta_X,delta_X)/np.dot(delta_X,delta_g))-(np.outer(np.dot(H,delta_g),delta_X)+np.outer(np.dot(H,delta_g),delta_X).transpose())/np.dot(delta_g,delta_X)
        H += dH

        err = LA.norm(g1)
        k += 1
        if((k%print_rate)==0):
            L2error = obj_func(X)
            print("k: %i\t err: %.3e\t L2error: %.6e" % (k,err,L2error))
        # store values
        appendToEachFile(err,L2error,k,X)
        # iter_store[k] = k
        # err_store[k] = err
        # L2error_store[k] = L2error
        # sol_store[k,:] = X

        # update for next iteration
        g0 = 1.0*g1

    kmax = k
    # sol_store = sol_store[:(kmax+1),:]
    # err_store = err_store[:(kmax+1)]
    # iter_store = iter_store[:(kmax+1)]
    # L2error_store = L2error_store[:(kmax+1)]
#=====================================================


#=====================================================
#                       MAIN
#=====================================================
#-----------------------------------------------------
global Nel, xL, xR, nVert
#-----------------------------------------------------
global L2error_file, vert_file, err_file, iter_file, data_dir
mdo_data_dir = "./Data/MDO/" # root directory for MDO data
#-----------------------------------------------------
global SetupFile, ParametersFile, VerticesFile, VerticesFile_sh
SetupFile="./SolverSetupFiles/MDO.setup"
ParametersFile = "./parameters/arctangent.parameters"
VerticesFile = "vertices.txt"
VerticesFile_sh = "./" + VerticesFile
#-----------------------------------------------------
global restart_switch, iter_restart
restart_level = 0 # USER INPUT; select 0 for not restarting
# Set the restart switch accordingly
if(restart_level==0):
    restart_switch = 0 # Not restarting
else:
    restart_switch = 1 # Yes restarting

if(restart_switch==0):
    #-----------------------------------------------------
    # FOR STARTING THE PROBLEM WITHOUT RESTART
    #-----------------------------------------------------
    print("---------------------------------------------")
    print("STARTING THE PROBLEM WITHOUT RESTART")
    print("---------------------------------------------")
    # Set the data directory
    data_dir = mdo_data_dir
    # Read setup file to obtain Nel, xL, xR
    readSetupFile(SetupFile)
    nVert = Nel-1
    # Initialize the vertices as equispaced
    xgrid = np.linspace(xL,xR,Nel+1)
    xvert_start = xgrid[1:Nel]
    writeVerticesToFile(xvert_start,VerticesFile) # NOT NEEDED
    
    # Set the initial vertices for optimization algorithm
    xvert = 1.0*xvert_start

    # Create the files for storing stuff
    L2error_file = open(mdo_data_dir+"L2error_MDO_store.txt",'w')
    vert_file = open(mdo_data_dir+"vertices_MDO_store.txt",'w')
    err_file = open(mdo_data_dir+"err_MDO_store.txt",'w')
    iter_file = open(mdo_data_dir+"iter_MDO_store.txt",'w')

    # Close the files for storing stuff
    L2error_file.close()
    err_file.close()
    vert_file.close()
    iter_file.close()

    plot_parameters_file = open(mdo_data_dir+"plot_parameters.txt",'w')
    appendIntToFile(nVert,plot_parameters_file)
    appendDoubleToFile(xL,plot_parameters_file)
    appendDoubleToFile(xR,plot_parameters_file)
    plot_parameters_file.close()
elif(restart_switch==1):
    #-----------------------------------------------------
    # FOR STARTING THE PROBLEM WITH RESTART
    #-----------------------------------------------------
    # Load important setup parameters: nVert, xL, xR
    nVert, xL, xR = np.loadtxt(mdo_data_dir+"plot_parameters.txt",unpack=True)
    nVert = int(nVert)

    # Restart directory
    restart_dir = mdo_data_dir
    for i in range(0,restart_level):
        restart_dir += "restart/"

    # Set the data directory
    data_dir = restart_dir

    # Load the restart point
    iter_restart = np.loadtxt(restart_dir+"iter_restart.txt")
    iter_restart = int(iter_restart)
    xvert_restart = np.loadtxt(restart_dir+"xvert_restart.txt")
    err_restart = np.loadtxt(restart_dir+"err_restart.txt")
    L2error_restart = np.loadtxt(restart_dir+"L2error_restart.txt")

    # Set the initial vertices for optimization algorithm
    xvert = 1.0*xvert_restart

    print("---------------------------------------------")
    print("RESTARTING COMPUTATION FROM:")
    print("- - - - - - - - - - - - - - - - - - - - - - -")
    print(" Restart level: %i" % restart_level)
    print(" Iteration: %i" % iter_restart)
    print(" Gradient norm: %.16e" % err_restart)
    print(" L2-error: %.16e" % L2error_restart)
    print("---------------------------------------------")
# Make Hessian storage directory
global H_store_directory
H_store_directory = data_dir + "hessian_store/"
subprocess.call(["mkdir",H_store_directory])
#-----------------------------------------------------

#-----------------------------------------------------
#               OPTIMIZATION SOLVER
#-----------------------------------------------------
global sol_store, err_store, iter_store, numvars, maxIter, print_rate, tol
numvars = 1*nVert
tol = 1.0e-6
maxIter = 100000
print_rate = 1

# Initial H and X0 for the optimization algorithm
X0 = 1.0*xvert
if(restart_switch==0):
    H = np.eye(numvars)
elif(restart_switch==1):
    H_restart_dir = mdo_data_dir
    for i in range(0,(restart_level-1)):
        H_restart_dir += "restart/"
    H_restart_dir += "hessian_store/"
    H_file_name = H_restart_dir + "H_iter_%i" % iter_restart
    H_file_name += ".txt"
    H = np.loadtxt(H_file_name)

# Optimization Algorithm
QuasiNewtonMethod(X0,H)
#=====================================================