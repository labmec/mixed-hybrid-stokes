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
#include "TPZIdentitySolver.h"

#include <TPZPardisoSolver.h>

const int global_nthread = 16;

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh, std::string filename, ProblemData *problem_data, TPZGeoMesh *gmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name);
void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh, ProblemData *problem_data);

int main()
{
#ifdef PZ_LOG
//    TPZLogger::InitializePZLOG("log4cxx.cfg");
    TPZLogger::InitializePZLOG();
#endif
    
    bool printdata = false;
    
    std::string filepath = "/Users/CarlosPuga/programming/HybridStokesResearch/DataInput/";
    std::string filename = "DifferentObstructions";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = simData.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::ENone;
    flow->fDimension = simData.Dim();
    
    // creating the meshes
    // geometric mesh
    TPZGeoMesh *gmesh = nullptr;
    gmesh = TPZMeshOperator::CreateGMesh(&simData, false);
    
    // velocity computational mesh
    TPZCompMesh *cmesh_v = nullptr;
    cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    // pressure computational mesh
    TPZCompMesh *cmesh_p;
    cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    // atomic meshes for static condensation
    TPZCompMesh *cmesh_Mp = nullptr; // mean pressure
    TPZCompMesh *cmesh_Mv = nullptr; // and mean flux
    
    if(simData.CondensedElements() && simData.HdivType() != ProblemData::EConstant)
    {
        cmesh_Mv = TPZMeshOperator::CreateCmeshG(&simData, gmesh);
        cmesh_Mp = TPZMeshOperator::CreateCmeshPm(&simData, gmesh);
    }

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = nullptr;
    cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh, flow);
    
    // static condensation
    if (simData.CondensedElements())
        TPZMeshOperator::CondenseElements(&simData, cmesh_m, gmesh);
    
    cmesh_m->SaddlePermute();
    
    RenumType renum = RenumType::EMetis;
    if (global_nthread == 0)
        renum = RenumType::ENone;
    
    TPZLinearAnalysis an(cmesh_m, renum);
    
    // starting the simulation
    SolveProblem(an, cmesh_m, filename, &simData, gmesh);
        
    // printing meshes
    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    // Post Process
    PrintResults(an, cmesh_m, &simData, filename);
    
    // Calculating error
    if (flow->fExactSol)
        EvaluateErrors(filename, an, flow, cmesh_m, &simData);
    
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

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, std::string filename, ProblemData *problem_data, TPZGeoMesh *gmesh)
{
    TPZSimpleTimer timer;
    
    TPZSSpStructMatrix<> strmat(cmesh_m);
    strmat.SetNumThreads(global_nthread);
        
    //    Applying filters to null normal and tangential velocities
    TPZEquationFilter filter(cmesh_m->NEquations());
    std::set<int64_t> removeEquations;

    //     on the boundary
    TPZMeshOperator::ConfigureBoundaryFilter(gmesh, cmesh_m, problem_data, removeEquations);

    //     between the domain and the obstruction
    if (problem_data->ObstructionID() != - 1)
        TPZMeshOperator::ConfigureObstructionFilter(gmesh, cmesh_m, problem_data, removeEquations);

    filter.ExcludeEquations(removeEquations);
    strmat.EquationFilter() = filter;
    
    an.SetStructuralMatrix(strmat);
    
    // setting solver
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
    
    TPZSimpleTimer solve_time("Solve time");
    
    an.Solve();
    
    // saving the simulation data (time of execution and number of equations)
    std::cout << "Time to solve = " << solve_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished solver...\n";
    
    std::cout << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << std::endl;
    
    simStatus << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << " s" << std::endl;
    simStatus << "Total NEquations: " << cmesh_m->Solution().Rows() << std::endl;
    simStatus << "Condensed NEquations: " << cmesh_m->NEquations() << std::endl;
    
    return;
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name)
{
    std::cout << "Starting PostProcess" << std::endl;
    TPZSimpleTimer postProc("Post proc time");
    
    // vtk export
    TPZVec<std::string> fields;
    
    if (problem_data->HasAnalyticSolution())
        fields =
        {
            "Pressure",
            "ExactPressure",
//            "ErrorPressure",
            
            "Velocity",
            "ExactVelocity",
//            "ErrorVelocity",
            
            "Stress",
            "ExactStress",
//            "ErrorStress"
        };
    else
        fields = {"Pressure", "Velocity", "Stress"};
    
    TPZVTKGenerator vtk(cmesh_m, fields, file_name, problem_data->Resolution());
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}

void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh_m, ProblemData *problem_data)
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
    an.SetThreadsForError(global_nthread);
    
    TPZMaterial *mat = cmesh_m->FindMaterial(problem_data->DomainVec()[0].matID);
    TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
    
    bool store_errors = false;
    std::ofstream ErrorOut(file_name  + "_Errors.txt");
    
    TPZSimpleTimer post_time("PostProcess time");
    
    an.PostProcessError(Errors, store_errors);
    
    std::cout << "Time PostProc Error = " << post_time.ReturnTimeDouble()/1000 << "s" << std::endl;
    
    ErrorOut << "###### Computed Errors ######" << std::endl;
    
    for (int i = 0; i < Errors.size(); i++)
        ErrorOut << "L2 " << fields[i] << " = " << Errors[i] << std::endl;
}
