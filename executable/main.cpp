#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include <TPZVTKGenerator.h>
#include <TPZSimpleTimer.h>

#include "ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"
#include <pzskylstrmatrix.h>
#include <pzfstrmatrix.h>
#include <pzstrmatrixot.h>
#include <TPZSpStructMatrix.h>

int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
//    TPZLogger::InitializePZLOG();
#endif
    bool printdata = false;
    int nThreads = 16;
    int nThreadsError = 16;
    
    std::string filepath = "/Users/CarlosPuga/programming/HybridStokesResearch/DataInput/";
    std::string filename = "Square_16_4_HdivS";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = simData.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::EPaperComp;
    flow->fDimension = simData.Dim();
    
    // creating the meshes
    // geometric mesh
    TPZGeoMesh *gmesh = TPZMeshOperator::CreateGMesh(&simData, false);
    
    // velocity computational mesh
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    // pressure computational mesh
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    // atomic meshes for static condensation
    TPZCompMesh *cmesh_Mp = nullptr; // mean pressure
    TPZCompMesh *cmesh_Mv = nullptr; // and mean flux
    if(simData.CondensedElements() && simData.HdivType() != ProblemData::EConstant)
    {
            cmesh_Mp = TPZMeshOperator::CreateCmeshPm(&simData, gmesh);
            cmesh_Mv = TPZMeshOperator::CreateCmeshG(&simData, gmesh);
    }

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh, flow);
    
    // static condensation
    if (simData.CondensedElements())
        TPZMeshOperator::CondenseElements(&simData, cmesh_m, gmesh);
    
    cmesh_m->SaddlePermute();
    
    TPZLinearAnalysis an(cmesh_m, RenumType::EMetis);
    
//    TPZSSpStructMatrix<STATE, TPZStructMatrixOT<STATE>> strmat(cmesh_m);
    TPZSSpStructMatrix<> strmat(cmesh_m);
//    TPZFStructMatrix<> strmat(cmesh_m);
//     TPZSkylineStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(nThreads);
    
//     Applying filters to null normal and tangential velocities
    TPZEquationFilter filter(cmesh_m->NEquations());
    std::set<int64_t> removeEquations;

//     on the boundary
    TPZMeshOperator::ConfigureBoundaryFilter(gmesh, cmesh_m, &simData, removeEquations);

//     between the domain and the obstruction
    if (simData.ObstructionID() != - 1)
        TPZMeshOperator::ConfigureObstructionFilter(gmesh, cmesh_m, &simData, removeEquations);

    filter.ExcludeEquations(removeEquations);
    strmat.EquationFilter() = filter;
    
    // starting the simulation
    TPZSimpleTimer timer;
    
    an.SetStructuralMatrix(strmat);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    std::ofstream simStatus(filename + "_Data.txt");
    
    // ... the assemble
    std::cout << "Starting assemble...\n";
    TPZSimpleTimer ass_time("Assemble time");
    an.Assemble();
    std::cout << "Time to assemble = " << ass_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished assemble...\n";
    
    // ... and the solver
    std::cout << "Starting solver...\n";
    std::cout << "Number of Equations " << cmesh_m->NEquations() << std::endl;
    
    auto matrix = an.MatrixSolver<STATE>().Matrix();
    
    TPZSimpleTimer solve_time("Solve time");
    an.Solve();
    std::cout << "Time to solve = " << solve_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished solver...\n";
    
    std::cout << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << std::endl;
    
    simStatus << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << " s" << std::endl;
    simStatus << "Total NEquations: " << cmesh_m->Solution().Rows() << std::endl;
    simStatus << "Condensed NEquations: " << cmesh_m->NEquations() << std::endl;
    
    // printing meshes
    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    // vtk export
    TPZVec<std::string> fields;
    
    if (flow->fExactSol)
        fields = 
        {
            "Pressure",
            "ExactPressure",
            "ErrorPressure",
            
            "Velocity",
            "ExactVelocity",
            "ErrorVelocity",
            
            "Stress",
            "ExactStress",
//            "ErrorStress"
        };
    else
     fields =
        {
            "Pressure",
        
            "Velocity",
        
            "Stress",
        };
    
    TPZVTKGenerator vtk(cmesh_m, fields, filename, 3);
    vtk.Do();
    
    // Calculating error
    
    if (flow->fExactSol != 0)
    {
        std::cout << "Calculating error..." << std::endl;
        TPZVec<std::string> fields =
        {
            "p error",
            "p_ex error",
            
            "u error",
            "u_ex error",
            
            "div error",
            "div_ex error",
            
            "sigma error",
            "sigma_ex error",
            
            "dev_sigma error",
            "dev_sigma_ex error"
        };
        
        an.SetExact(flow->ExactSolution());
        an.SetThreadsForError(nThreadsError);
        
        TPZMaterial *mat = cmesh_m->FindMaterial(simData.DomainVec()[0].matID);
        TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
        TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
        
        bool store_errors = false;
        std::ofstream ErrorOut(filename  + "_Errors.txt");
        
        TPZSimpleTimer post_time("PostProcess time");
        an.PostProcessError(Errors, store_errors);
        std::cout << "Time PostProc Error = " << post_time.ReturnTimeDouble()/1000 << "s" << std::endl;
        
        ErrorOut << "###### Computed Errors ######" << std::endl;
        for (int i = 0; i < Errors.size(); i++)
            ErrorOut << "L2 " << fields[i] << " = " << Errors[i] << std::endl;
    }
    
    // Deleting pointers
    if (cmesh_m)
        delete cmesh_m;
    
    if (cmesh_v)
        delete cmesh_v;
    
    if (cmesh_p)
        delete cmesh_p;
    
    if (cmesh_Mp)
        delete cmesh_Mp;
    
    if (cmesh_Mv)
        delete cmesh_Mv;
    
    if (gmesh)
        delete gmesh;
    
    std::cout << "\n\nSimulation finished without errors :) \n\n";
    
	return 0;
}
