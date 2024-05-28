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

// **********************
// FUNCTION DECLARATIONS
// **********************
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, std::string filename);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name);
void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh_m, ProblemData *problem_data);

const int global_nThreads = 16;

// **********************
//     MAIN FUNCTION
// **********************
int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    
    bool print_mesh = false; // turn this on to print the mesh files
    bool blend = true; // switch this to use or not blend elements
    
    std::string blend_file = (blend) ? "_B" : "_NB";
    
    // reading problem data from json file
    std::string file_path = "/Users/CarlosPuga/programming/HybridStokesResearch/DataInput/"; // CHANGE THIS
    std::string file_name = "TaylorCouette_1_2_S"; // first number is the VelpOrder, second number is the refinement level. Letter indicates whether it is Hdiv Standard or Constant
    
    ProblemData problem_data;
    
    problem_data.ReadJson(file_path + file_name + ".json");
    
    // creating an analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = problem_data.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fDimension = problem_data.Dim();
    flow->fExactSol = TStokesAnalytic::ETaylorCouette;
    
    // creating gmsh
    TPZGeoMesh *gmesh = TPZMeshOperator::CreateGMesh(&problem_data, blend);
        
    // creating velocity cmesh
    TPZCompMesh *cmesh_v = TPZMeshOperator::CreateCMeshV(&problem_data, gmesh);
    
    // creatiing pressure cmesh
    TPZCompMesh *cmesh_p = TPZMeshOperator::CreateCmeshP(&problem_data, gmesh);
    
    // atomic meshes for static condensation
    TPZCompMesh *cmesh_Pm = nullptr; // mean pressure
    TPZCompMesh *cmesh_G = nullptr; // and mean flux
    
    if (problem_data.CondensedElements() && problem_data.HdivType() != ProblemData::EConstant)
    {
        cmesh_G = TPZMeshOperator::CreateCmeshG(&problem_data, gmesh);
        cmesh_Pm = TPZMeshOperator::CreateCmeshPm(&problem_data, gmesh);
    }
    
    // multiphysics cmesh
    TPZMultiphysicsCompMesh *cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&problem_data, gmesh, flow);
    
    // static condensation
    if (problem_data.CondensedElements())
        TPZMeshOperator::CondenseElements(&problem_data, cmesh_m, gmesh);
    
    cmesh_m->SaddlePermute();
    
    // Analysis
    // Solve Multiphysics
    RenumType renum = RenumType::EMetis;
    if (global_nThreads == 0)
        renum = RenumType::ENone;
        
    TPZLinearAnalysis an(cmesh_m, renum);
    SolveProblemDirect(an, cmesh_m, file_name);
    
    // Post Process
    PrintResults(an, cmesh_m, &problem_data, file_name + blend_file);
    
    // Error
    if (flow->fExactSol)
        EvaluateErrors(file_name + blend_file, an, flow, cmesh_m, &problem_data);
    
    // Deleting mesh pointers
    if (cmesh_m)
        delete cmesh_m;

    if (cmesh_v)
        delete cmesh_v;

    if (cmesh_p)
        delete cmesh_p;
    
    if (cmesh_Pm)
        delete cmesh_Pm;
    
    if (cmesh_G)
        delete cmesh_G;
    
    if (gmesh)
        delete gmesh;
    
    std::cout << "--------- Simulation finished ---------" << std::endl;
    
    return 0;
    
}

// **********************
//    SOLVER FUNCTION
// **********************
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, std::string filename)
{
    TPZSimpleTimer timer;
    
    TPZSSpStructMatrix<> strmat(cmesh_m);
    strmat.SetNumThreads(global_nThreads);
    an.SetStructuralMatrix(strmat);
    
    // seting solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    // assembling the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_assemble;
    an.Assemble();
    std::cout << "Time to assemble = " << time_assemble.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    // solving the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_solver;
    an.Solve();
    std::cout << "Time to solve = " << time_solver.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    std::cout << "Simulation time = " << timer.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}

// **********************
// POST PROCESS FUNCTION
// **********************
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name)
{
    std::cout << "--------- PostProcess ---------" << std::endl;
    TPZSimpleTimer time_postProc;
    
    TPZVec<std::string> fields = {"Pressure", "Velocity", "Stress"};
    
    auto vtk = TPZVTKGenerator(cmesh_m, fields, file_name, problem_data->Resolution());
    vtk.SetNThreads(global_nThreads);
    vtk.Do();
    
    std::cout << "Total time = " << time_postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}

void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh_m, ProblemData *problem_data)
{
    TPZSimpleTimer time_error;
    
    std::cout << "--------- Error  ---------" << std::endl;
    TPZVec<std::string> fields = {
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
    an.SetThreadsForError(global_nThreads);
    
    TPZMaterial *mat = cmesh_m->FindMaterial(problem_data->DomainVec()[0].matID);
    TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(mat);
    TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
    
    bool store_errors = false;
    std::ofstream ErrorOut(file_name + "_Errors.txt");
    
    an.PostProcessError(Errors, store_errors);
    std::cout << "Error time = " << time_error.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    ErrorOut << "###### Computed Errors ######" << std::endl;
    for (int i = 0; i < Errors.size(); i++)
        ErrorOut << "L2 " << fields[i] << " = " << Errors[i] << std::endl;
    
    return;
}
