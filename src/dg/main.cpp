#include <iostream>
#include <fstream> // for writing to files
#include <string> // for strings
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <sstream> // for user input from command line
// polynomial library from NEKTAR
#include "polylib.h"
// My standard libraires:
#include "allocatePointersLib.h"
#include "writeToFile.h"
#include "flattenLib.h"
#include "strManipLib.h"
// My Spectral Element Method (SEM) libraries
#include "var.h"
#include "allocateVar.h"
#include "dglib.h"
#include "ErrorEstimates.h"
#include "PDE.h"
#include "readFilesDG.h"
#include "initZeros.h"

using namespace polylib;
using namespace std;

//===============================================
//              LAPACK FUNCTIONS
//===============================================
extern "C" {extern void dgesv_(int *, int *, double (*), int *, int [], double (*), int *, int*);}
//===============================================
//  DEFINE FUNCTIONS not classified to a .h file
//===============================================
int getNDUMP(int NOUT, int DUMP_INTERVAL)
{
    int i_dump = 0;
    for(int i=0; i<NOUT; i++)
    {
        if((i % DUMP_INTERVAL)==0)
        {
            i_dump += 1;
        }
    }
    return i_dump;
}
//***********************************************
//***********************************************
//***********************************************
int main()
{
    // User input Nel and P
    string cmdline_userinput; // stores what is entered in cmd line
    cout << endl;
    cout << "USER INPUT:" << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "Convergence study?\n\t(0) No\n\t(1) Yes\n";
    getline(cin,cmdline_userinput);
    stringstream(cmdline_userinput) >> ConvStudyFlag;
    if(ConvStudyFlag == 1)
    {
        cout << "Enter polynomial order: P = ";
        getline(cin,cmdline_userinput);
        stringstream(cmdline_userinput) >> P;
        cout << "Enter number of elements: Nel = ";
        getline(cin,cmdline_userinput);
        stringstream(cmdline_userinput) >> Nel;
        cout << "Select type of study:\n\t(0) p-convergence\n\t(1) h-convergence\n";
        getline(cin,cmdline_userinput);
        stringstream(cmdline_userinput) >> ConvStudyType;
    }
    else if(ConvStudyFlag == 0)
    {
        cout << "Store transient?\n\t(0) No\n\t(1) Yes\n";
        getline(cin,cmdline_userinput);
        stringstream(cmdline_userinput) >> store_transient;
    }
    cout << "-----------------------------------------------" << endl;

    //-----------------------------------------------
    /* Read the setup file name */
    getline(cin,cmdline_userinput);
    namesetup = cmdline_userinput;
    readSetupFile(namesetup);
    //-----------------------------------------------

    if(problem=="SteadyStateLinearAdvection_Manufactured")
    {
        LinAdvSpeed = 1.0;
    }
    else
    {
        LinAdvSpeed = 2.0*PI; // for LinearAdvection, irrelevant for Burgers
    }
    /* Compute number of points, problem size, etc */
    updateNumPoints_ProblemSizeParameters(P, Nel);
    /* Display the parameters */
    displayMainParameters();
    /* Allocate all variables (except xDomainBounds) */
    allocateVariables();
    
    //===============================================
    //      Integeration & Differentiation
    //===============================================
    //-----------------------------------------------
    //          Gauss-Lobatto-Legendre
    //-----------------------------------------------
    zwgll(z, w, np); // get zeros and weights for GLL quadrature
    Dgll(D_flat, Dt_flat, z, np); // get differentiation matrix for GLL quadrature (C-flattened)
    unflattenC(np, np, D, D_flat); // unflatten D matrix
    //===============================================
    //          One time calculations / Inits
    //===============================================
    //-----------------------------------------------
    //                  Grid
    //-----------------------------------------------
    if(subproblem=="MDO")
    {
        /* Read the setup file name */
        getline(cin,cmdline_userinput);
        // namevertices = cmdline_userinput;
        readElementVertices(cmdline_userinput);
        //-----------------------------------------------
        cout << "ELEMENT BOUNDS -- Python" << endl;
        for(int i=0; i<(Nel+1); i++)
        {
            cout << elBounds[i] << endl;
        }
    }
    else
    {
        // if not doing MDO
        elementBounds(Nel, xDomainBounds, elBounds);
        // cout << "ELEMENT BOUNDS" << endl;
        // for(int i=0; i<(Nel+1); i++)
        // {
        //  cout << elBounds[i] << endl;
        // }
    }
    //-----------------------------------------------
    if(problem=="SteadyStateLinearAdvection_Manufactured")
    {
        /* Read the setup file name */
        getline(cin,cmdline_userinput);
        // nameArctanParam = cmdline_userinput;
        readArctanParameters(cmdline_userinput);
        cout << "arctan_const: " << arctan_const << endl;
        cout << "arctan_magnitude: " << arctan_magnitude << endl;
        cout << "arctan_NumShocks: " << arctan_NumShocks << endl;
        for(int i=0; i<arctan_NumShocks; i++)
        {
            cout << "arctan_ShockMag[" << i << "]: " << arctan_ShockMag[i] << endl;
            cout << "arctan_ShockLoc[" << i << "]: " << arctan_ShockLoc[i] << endl;
        }
    }

    // return 0; // for checking the reading and vertex input

    /* One time calculation -- setup */
    for(int eln=0; eln<Nel; eln++)
    {
        /* Setup grid */
        chi(np, eln, z, elBounds, x, Jac); // map to element (x == z_mapped) -- i.e. x=chi(z) where z=xi
        xGlob(np, Nel, eln, x, xG);
    }
    /* Store elementwise grid */
    writeMatrixToFile(Nel,np,x,"x","txt","Data");
    //-----------------------------------------------
    //                  Mapping
    //-----------------------------------------------
    mapping_frequency(Nel, P, map); // delete later on // useless?
    mapping_frequency(Nel, P, map_frequency); // used for global mapping
    mapping_physical(Nel, np, map_physical);
    //-----------------------------------------------
    //          One time initializations
    //-----------------------------------------------   
    basis_vec_xLR(np, P, z, phi1, phi_L, phi_R);
    for(int eln=0; eln<Nel; eln++)
    {
        /* Compute M,C,L matrices for each element */
        assem(Mdim, np, P, eln, z, w, phi1, phi2, dphi1, dphi2, Jac, D, M, C, L); // assemble mass matrix
    }

    // for(int eln=0; eln<Nel; eln++)
    // {
    //  cout << "eln: " << eln << "\t xL = " << xG[map_physical[eln][0]] << "\t xR = " << xG[map_physical[eln][np-1]] << endl;
    // }

    if(problem!="SteadyStateLinearAdvection_Manufactured")
    {
        /* Initial condition */
        for(int eln=0; eln<Nel; eln++)
        {
            /* Initial condition on physical grid */
            InitialCondition(np, eln, x, u0);
            /* Initial condition in the frequency space */
            PhysicalToFrequency(P, eln, Jac, u0, u_hat);
            // for(int i=0; i<np; i++)
            // {
            //  cout << "eln: " << eln << "\t i: " << i << "\t u0: " << u0[eln][i] << endl;
            // }
        }
        
        /* Double check initialization of the solution with IC */
        for(int eln=0; eln<Nel; eln++)
        {
            FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_delta);
        }

        /* Store initialization in transient format */
        if(store_transient == 1)
        {
            iout = 1;
            writeTransientSolutionToFile(Nel, np, u_delta, "u_delta", "txt", "Data/Transient", iout, T0);
            writeTransientSolutionToFile(Nel, np, u0, "u_exact", "txt", "Data/Transient", iout, T0);
        }
        else if(store_transient == 0)
        {
            writeMatrixToFile(Nel,np,u0,"u0","txt","Data");
            writeMatrixToFile(Nel,np,u_delta,"u0_delta","txt","Data");
        }

        /* Compute smallest delta x */
        DELTA_X_MIN = xG[1]-xG[0]; // initialize
        double dx;
        for(int i=2; i<npG; i++)
        {
            dx = xG[i]-xG[i-1];
            if(dx < DELTA_X_MIN)
            {
                DELTA_X_MIN = dx;
            }
        }
        // writeArrayToFile(np,z,"z","txt","Data");
        // writeArrayToFile(npG,xG,"xG","txt","Data");
        cout << "DELTA_X_MIN: " << DELTA_X_MIN << endl;
        /* Compute time step based on CFL number */
        if(problem=="LinearAdvection")
        {
            TIME_STEP = (CFL*DELTA_X_MIN)/LinAdvSpeed;
        }
        else if((problem=="NonhomogenousInviscidBurgers") || (problem=="HomogenousInviscidBurgers"))
        {
            TIME_STEP = getTimeStep_Burgers(u_hat);
        }

        // double *t_store = dvector(NDUMP);

        double t = T0;
        int NOUT = 0;

        while(t<Tf)
        {
            t += TIME_STEP;
            NOUT += 1;
        }
        cout << "CFL: " << CFL << endl;
        cout << "NOUT: " << NOUT << endl;
        cout << "t: " << t << endl;
        cout << "TIME_STEP: " << TIME_STEP << endl;

        /* Correct time step for exactly reaching Tf */
        // - This will be slightly smaller
        if(problem=="LinearAdvection")
        {   
            TIME_STEP = (Tf-T0)/double(NOUT);
            CFL = LinAdvSpeed*TIME_STEP/DELTA_X_MIN; // update CFL
        }

        t = T0;
        NOUT = 0;
        while(abs(t-Tf)>(1.0e-10)) // while(t<Tf)
        {
            if((problem == "NonhomogenousInviscidBurgers") || (problem == "HomogenousInviscidBurgers"))
            {
                TIME_STEP = getTimeStep_Burgers(u_hat);
                
                if((t+TIME_STEP)>Tf)
                {
                    TIME_STEP = Tf - t;
                }
            }
            NOUT += 1;
            RK4(t, TIME_STEP, u_hat);
            t += TIME_STEP;
            // cout << "iout: " << NOUT << "\t" << "t: " << t << "\t TIME_STEP: " << TIME_STEP << endl;

            /* Storing the transient based on dump interval */
            if(store_transient == 1)
            {
                if((NOUT % DUMP_INTERVAL)==0)
                {
                    iout += 1;
                    /* Transfer to the phyiscal space */
                    for(int eln=0; eln<Nel; eln++)
                    {
                        FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_delta);
                        ExactSolution(np, eln, t, x, u_exact);
                    }
                    writeTransientSolutionToFile(Nel, np, u_delta, "u_delta", "txt", "Data/Transient", iout, t);
                    writeTransientSolutionToFile(Nel, np, u_exact, "u_exact", "txt", "Data/Transient", iout, t);
                }
            }
        }

        cout << "CFL: " << CFL << endl;
        cout << "NOUT: " << NOUT << endl;
        cout << "t: " << t << endl;
        cout << "TIME_STEP: " << TIME_STEP << endl;

        /* Final solution: Transfer to the phyiscal space */
        for(int eln=0; eln<Nel; eln++)
        {
            FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_delta);
            ExactSolution(np, eln, t, x, u_exact);
        }

        /* Store final solution in transient format */
        if(store_transient == 1)
        {
            iout += 1;
            writeTransientSolutionToFile(Nel, np, u_delta, "u_delta", "txt", "Data/Transient", iout, t);
            writeTransientSolutionToFile(Nel, np, u_exact, "u_exact", "txt", "Data/Transient", iout, t);
            writeIntToFile(iout, "NDUMP","txt", "Data/Transient");
        }
        else if(store_transient == 0)
        {
            writeMatrixToFile(Nel,np,u_exact,"u_exact","txt","Data");
            writeMatrixToFile(Nel,np,u_delta,"uF_delta","txt","Data");
        }
    }
    else
    {
        cout << "-----------------------------------------------" << endl;
        cout << "\t STEADY STATE IMPLICIT " << endl;
        cout << "-----------------------------------------------" << endl;
            
        for(int eln=0; eln<Nel; eln++)
        {
            /* Compute source term */
            SourceTerm(np, eln, 0.0, x, q_source);
            /* Compute inner product of source term and phi */
            InnerProduct_f_phi(np, P, eln, phi1, Jac, q_source, q_phi_innerproduct);
        }

        /* Initialize as zeros */
        initZeros_dvector(GlobalDim, u_hat_global);
        initZeros_dvector(GlobalDim, RHS_VEC_global);
        initZeros_dvector(GlobalDim, q_phi_innerproduct_global);
        initZeros_dvector(GlobalDim*GlobalDim, LHS_MAT_global_FLAT);
        initZeros_dmatrix(GlobalDim, GlobalDim, LHS_MAT_global);
        initZeros_dmatrix(GlobalDim, GlobalDim, M_global);
        initZeros_dmatrix(GlobalDim, GlobalDim, C_global);
        initZeros_dmatrix(GlobalDim, GlobalDim, L_global);
        initZeros_dmatrix(Mdim, Mdim, PHI_LL_MAT);
        initZeros_dmatrix(Mdim, Mdim, PHI_LR_MAT);
        initZeros_dmatrix(GlobalDim, GlobalDim, PHI_MAT_global);

        /* Transform the elementwise standard matrices to global */
        ElementwiseToGlobal_Matrix(P, Nel, map_frequency, M, M_global);
        ElementwiseToGlobal_Matrix(P, Nel, map_frequency, C, C_global);
        ElementwiseToGlobal_Matrix(P, Nel, map_frequency, L, L_global);
        /* Transform the elementwise q_phi_innerproduct (source term) to global */
        ElementwiseToGlobal_Vector(P, Nel, map_frequency, q_phi_innerproduct, q_phi_innerproduct_global);
        /* Compute the PHI_LL and PHI_LR square matrices of dimension P+1 */
        compute_PHI_MATS();
        /*  Assemble global PHI matrix (strong form) */
        assemGlobalPhiMatrix();
        /* Assemble the linear system to be solved */
        assemGlobalMatrixLHS(LHS_MAT_global);
        assemGlobalVectorRHS(RHS_VEC_global);

        double **C_1el = dmatrix(Mdim,Mdim);
        for(int i=0; i<Mdim; i++)
        {
            for(int j=0; j<Mdim; j++)
            {
                C_1el[i][j] = C[i][j][0];
            }
        }

        /* Write to files for checking matrices and DG setup */
        // writeMatrixToFile(Nel,Mdim,q_phi_innerproduct,"q_phi_innerproduct","txt","Data");
        // writeArrayToFile(GlobalDim,q_phi_innerproduct_global,"q_phi_innerproduct_global","txt","Data");
        // writeMatrixToFile(Mdim,Mdim,PHI_LL_MAT,"PHI_LL_MAT","txt","Data");
        // writeMatrixToFile(Mdim,Mdim,PHI_LR_MAT,"PHI_LR_MAT","txt","Data");
        // writeMatrixToFile(GlobalDim,GlobalDim,PHI_MAT_global,"PHI_MAT_global","txt","Data");
        // writeMatrixToFile(Mdim,Mdim,C_1el,"C_1el","txt","Data");
        // writeMatrixToFile(GlobalDim,GlobalDim,C_global,"C_global","txt","Data");
        // writeMatrixToFile(GlobalDim,GlobalDim,LHS_MAT_global,"LHS_MAT_global","txt","Data");
        // writeArrayToFile(GlobalDim,RHS_VEC_global,"RHS_VEC_global","txt","Data");

        /* Solve linear system to obtain the weights */
        /* - (1) Flatten matrix for LAPACK */
        flattenF(GlobalDim, GlobalDim, LHS_MAT_global, LHS_MAT_global_FLAT); // no issues here
        /* - (2) Solve for weights using LAPACK --> Stored in 6th argument */
        dgesv_(&GlobalDim, &NRHS, LHS_MAT_global_FLAT, &GlobalDim, ipiv_global, RHS_VEC_global, &GlobalDim, &INFO);

        /* - (3) Store the solution (weights) for this element */
        for(int i=0; i<GlobalDim; i++)
        {
            u_hat_global[i] = RHS_VEC_global[i];
        }
        /* Tranform global solution vector to elementise vector */
        GlobalToElementwise_Vector(P, Nel, map_frequency, u_hat_global, u_hat);
        /* Output u_hat */
        // writeMatrixToFile(Nel,Mdim,u_hat,"u_hat","txt","Data");
        // writeArrayToFile(GlobalDim,u_hat_global,"u_hat_global","txt","Data");
        /* Transfer to the phyiscal space + Compute exact solution */
        for(int eln=0; eln<Nel; eln++)
        {
            FrequencyToPhysical(np, P, eln, z, phi1, u_hat, u_delta);
            ExactSolution(np, eln, 0.0, x, u_exact);
        }
        writeMatrixToFile(Nel,np,u_exact,"u_exact","txt","Data");
        writeMatrixToFile(Nel,np,u_delta,"uF_delta","txt","Data");
        cout << "LinAdvSpeed: " << LinAdvSpeed << endl;
    }

    double L2error = 0.0;
    L2error = L2_NORM(u_delta,u_exact);
    cout << "L2error: " << L2error << endl;
    writeDoubleToFile(L2error,"L2error","txt","Data");

    /* Append error estimate to appropriate file */
    std::ofstream logging;
    string fileName_ConvStudy;
    if(ConvStudyFlag == 1)
    {
        if(ConvStudyType == 0)
        {
            fileName_ConvStudy = "p_type_extension.txt";
            
        }
        else if(ConvStudyType == 1)
        {
            fileName_ConvStudy = "h_type_extension.txt";
        }

        logging.open(fileName_ConvStudy, std::ios_base::app);
        logging << nDOF << "\t" << L2error << "\t" << CFL << "\t" << P << "\t" << Nel << "\t" << np <<"\n";
        logging.close();
    }
    

    // Delete all pointers ?
    // delete [] z;delete [] w;delete [] p;delete [] elBounds;delete [] x;delete [] Jac;
    // delete [] phi1;delete [] phi2;delete [] f;delete [] u_delta;delete [] ipiv;
    // delete [] D;delete [] Dt;delete [] M;delete [] MG;delete [] fG;delete [] MG_flat;
    // delete [] xG;delete [] pG;delete [] D_flat;delete [] Dt_flat;delete [] dphi1;
    // delete [] dphi2;delete [] C;delete [] CG;delete [] CG_flat;delete [] L;delete [] LG;
    // delete [] LG_flat;delete [] gN;delete [] AG;delete [] gD;delete [] AG_flat;
    // delete [] fG_BC;delete [] xDomainBounds;
    
    return 0;
}