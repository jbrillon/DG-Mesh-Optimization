import numpy as np
import subprocess

def checkVertices(xVert):
    global xL, xR, nVert
    for i in range(1,nVert):
        if((xVert[i]>xVert[i-1]) or ((xVert[i]>xR) or (xVert[i]<xL))):
            print("WARNING: Overlapping vertice at x = %.5f" % xVert[i])

#=====================================================
#           SETUP THE RESTART FILES
#-----------------------------------------------------
def setup_restart_files(MDO_data_root_dir,restart_level):

    dir_to_restart_from = MDO_data_root_dir

    # Directory where the data currently sits
    for i in range(0,(restart_level-1)):
        dir_to_restart_from += "restart/"

    # Directory where the new restart files will be written
    dir_to_write_new_restart_files = dir_to_restart_from + "restart/"

    # Create the directory if it doesn't already exist
    subprocess.call(["mkdir", dir_to_write_new_restart_files])

    L2error_store = np.loadtxt(dir_to_restart_from+"L2error_MDO_store.txt")
    sol_store = np.loadtxt(dir_to_restart_from+"vertices_MDO_store.txt")
    err_store = np.loadtxt(dir_to_restart_from+"err_MDO_store.txt")
    iter_store = np.loadtxt(dir_to_restart_from+"iter_MDO_store.txt")

    # Find smallest optimization error (i.e. gradient norm) 
    # - (ignore last iter where it crashed since NaN)
    restart_index = np.argmin(abs(err_store[:-1]))-3
    fileline_number = restart_index+1

    # Trim the data up to restart index
    iter_store = iter_store[:(restart_index+1)]
    err_store = err_store[:(restart_index+1)]
    L2error_store = L2error_store[:(restart_index+1)]
    sol_store = sol_store[:(restart_index+1),:]

    # Create/open the new restart files
    L2error_file = open(dir_to_write_new_restart_files+"L2error_MDO_store.txt",'w')
    vert_file = open(dir_to_write_new_restart_files+"vertices_MDO_store.txt",'w')
    err_file = open(dir_to_write_new_restart_files+"err_MDO_store.txt",'w')
    iter_file = open(dir_to_write_new_restart_files+"iter_MDO_store.txt",'w')

    # Write the new restart files
    for i in range(0,restart_index+1):
        # Write iterations
        wstr = "%i\n" % int(iter_store[i])
        iter_file.write(wstr)

        # Write err
        wstr = "%.16e\n" % err_store[i]
        err_file.write(wstr)

        # Write L2error
        wstr = "%.16e\n" % L2error_store[i]
        L2error_file.write(wstr)

        # Write vertices to the file
        for v in (sol_store[i,:]):
            wstr = "%.16e" % v
            vert_file.write(wstr)
            wstr = " "
            vert_file.write(wstr)
        wstr = "\n"
        vert_file.write(wstr)
    
    # Close the new restart files
    vert_file.close()
    L2error_file.close()
    err_file.close()
    iter_file.close()

    # Files for easy restart
    restart_vert_file = open(dir_to_write_new_restart_files+"xvert_restart.txt",'w')
    restart_iter_file = open(dir_to_write_new_restart_files+"iter_restart.txt",'w')
    restart_err_file = open(dir_to_write_new_restart_files+"err_restart.txt",'w')
    restart_L2error_file = open(dir_to_write_new_restart_files+"L2error_restart.txt",'w')

    # Write iteration to start from
    wstr = "%i\n" % iter_store[restart_index]
    restart_iter_file.write(wstr)

    # Write err to start from
    wstr = "%.16e\n" % err_store[restart_index]
    restart_err_file.write(wstr)

    # Write L2error to start from
    wstr = "%.16e\n" % L2error_store[restart_index]
    restart_L2error_file.write(wstr)
    
    # Write vertices to start from
    for v in (sol_store[restart_index,:]):
        wstr = "%.16e" % v
        restart_vert_file.write(wstr)
        wstr = " "
        restart_vert_file.write(wstr)
    wstr = "\n"
    restart_vert_file.write(wstr)

    # Close the easy restart files
    restart_vert_file.close()
    restart_iter_file.close()
    restart_err_file.close()
    restart_L2error_file.close()
#=====================================================
restart_level = 1 # USER INPUT (could move it to a bash script)

if(restart_level>0):
    setup_restart_files("./Data/MDO/",restart_level)


# To be moved to the DG_MDO.py code
