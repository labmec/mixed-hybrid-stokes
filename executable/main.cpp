#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>
#include <TPZSSpStructMatrix.h>
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
    bool printdata = true;
    int nThreads = 16;
    
    std::string filepath = "../DataInput/";
    std::string filename = "TaylorCouette1";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = simData.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::ETaylorCouette;
    
    // creating the meshes
    // geometric mesh
    TPZGeoMesh *gmesh = TPZMeshOperator::CreateGMesh(&simData);
    
    // velocity computational mesh
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    // pressure computational mesh
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    // atomic meshes for static condensation
    TPZCompMesh *cmesh_Mp = nullptr; // mean pressure
    TPZCompMesh *cmesh_Mv = nullptr; // and mean flux
    if(simData.CondensedElements() && simData.HdivType() != ProblemData::EConstant){
            cmesh_Mp = TPZMeshOperator::CreateCmeshPm(&simData, gmesh);
            cmesh_Mv = TPZMeshOperator::CreateCmeshG(&simData, gmesh);
    }

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh, flow);
    
    // static condensation
    if (simData.CondensedElements())
        TPZMeshOperator::CondenseElements(&simData, cmesh_m, gmesh);
    
    cmesh_m->SaddlePermute();
    
    TPZLinearAnalysis an(cmesh_m, RenumType::ENone);
    
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
    
    // printing meshes
    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    // starting the simulation
    TPZSimpleTimer timer;
    
    an.SetStructuralMatrix(strmat);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    // ... the assemble
    std::cout << "Starting assemble...\n";
    TPZSimpleTimer ass_time("Assemble time");
    an.Assemble();
    std::cout << "Time to assemble = " << ass_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished assemble...\n";
    
    // ... and the solver
    std::cout << "Starting solver...\n";
    std::cout << "Number of Equations " << cmesh_m->NEquations() << std::endl;
    
    TPZSimpleTimer solve_time("Solve time");
    an.Solve();
    std::cout << "Time to solve = " << solve_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished solver...\n";
    
    std::cout << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << std::endl;
    
    // printing meshes
    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    // vtk export
    TPZVec<std::string> fields = 
    {
        "Pressure",
        "PressExact",
        "PressElError",
        
        "Velocity",
        "VelExact",
        "VelElError",
        
//        "Stress",
//        "StressExact",
//        "StressElError"
    };
    
    TPZVTKGenerator vtk(cmesh_m, fields, filename, 2);
    vtk.Do();
    
    // Calculating error
    std::cout << "Calculating error..." << std::endl;
    if (flow->fExactSol != 0)
    {
        std::vector<std::string> fields = {"p error", "p_ex error", "u error", "u_ex error", "div error", "div_ex error", "sigma error", "sigma_ex error"};
        
        an.SetExact(flow->ExactSolution());
        an.SetThreadsForError(nThreads);
        
        TPZMaterial *mat = cmesh_m->FindMaterial(simData.DomainVec()[0].matID);
        TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
        TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
        
        bool store_errors = false;
        std::ofstream ErrorOut(filename+"_Errors.txt");
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        an.PostProcessError(Errors, store_errors);
        
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time PostProc Error = " << (end - begin).count()/1000. << "s" << std::endl;
        
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


