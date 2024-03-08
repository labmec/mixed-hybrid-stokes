#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "TPZNullMaterial.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZSimpleTimer.h"
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include "TPZGmshReader.h"
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzintel.h"
#include "TPZNullMaterialCS.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzskylmat.h"
#include "TPZMatrixSolver.h"

#include "TPZStokesMaterialTH.h"

#include "ProblemData.h"

const int global_nthread = 0;

// **********************
// FUNCTION DECLARATIONS
// **********************

TPZGeoMesh *ReadMeshFromGmsh(std::string file_name, ProblemData *problem_data);
TPZCompMesh *CreateCMeshV(ProblemData *problem_data, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshP(ProblemData *problem_data, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *CreateMultiphysicsCMesh(ProblemData *problem_data, TPZGeoMesh *gmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data, std::string file_name);
void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh, ProblemData *problem_data)

#ifdef PZ_LOG
static TPZLogger logger("pz.TaylorHood");
#endif

// **********************
//     MAIN FUNCTION
// **********************
int main(int argc, char *argv())
{
    std::cout << "--------- Starting simulation ---------" << std::endl;
    
    // reading problem data from json file
    std::string file = "CouetteFlowTH";
    
    if(argc > 1)
        file = std::string(argv[1]);
    
    ProblemData problem_data;
    problem_data.ReadJson(file);
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = problem_data.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::EPaperComp;
    flow->fDimension = problem_data.Dim();
    
    // create mesh
    TPZGeoMesh *gmesh = nullptr;
    std::string mesh_file = problem_data.MeshName();
    
    gmesh = ReadMeshFromGmsh(mesh_file, &problem_data);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
    }
    
    // create computational meshes
    problem_data.MeshVector().resize(2);
    
    if (problem_data.DomainVec().size() > 1)
        DebugStop(); // please, implement the next lines correctly if many domains
    
    // velocity
    TPZCompMesh cmesh_v = CreateCMeshV(&problem_data, gmesh);
    {
        std::ofstream out("cmesh_u.txt");
        cmesh_v->Print(out);
    }
    
    // pressure
    TPZCompMesh cmesh_p = CreateCMeshP(&problem_data, gmesh);
    {
        std::ofstream out("cmesh_p.txt");
        cmesh_p->Print(out);
    }
    
    // multiphysics
    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiphysicsCMesh(&problem_data, gmesh);
    
    std::cout << cmesh_m->NEquations() << std::endl;
    
    // Analysis
    // Solve Multiphysics
    RenumType renum = RenumType::EMetis;
    if (global_nthread == 0)
        renum = RenumType::ENone;
    
    TPZLinearAnalysis an(cmesh_m, renum);
    SolveProblemDirect(an, cmesh_m);
    {
        std::ofstream out("cmesh_m.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out);
        std::ofstream out2("cmesh_m.txt");
        cmesh_m->Print(out2);
    }
    
    // Post Process
    PrintResults(an, cmesh_m, &problem_data, file);
    
    // calculating error
    if (flow->fExactSol != 0)
        EvaluateErrors(file, an, flow, cmesh, &problem_data);
    
    // Deleting pointers
    if (cmesh_m)
        delete cmesh_m;
    
    if (cmesh_v)
        delete cmesh_v;
    
    if (cmesh_p)
        delete cmesh_p;
    
    if (gmesh)
        delete gmesh;
    
    std::cout << "--------- Simulation finished ---------" << std::endl;
}
// *        END          *


// **********************
//       GEO MESH
// **********************
TPZGeoMesh *ReadMeshFromGmsh(std::string file_name, ProblemData *problem_data)
{
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    
    TPZGmshReader reader;
    reader.GeometricGmshMesh(file_name, gmesh);
    
    gmesh->BuildConnectivity();
    
    return gmesh;
}
// *        END          *

// **********************
//    VEL COMP MESH
// **********************
TPZCompMesh *CreateCMeshV(ProblemData *problem_data, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v->SetName("CMesh_U");
    
    const int dimension = problem_data->Dim();
    std::set<int> materialIDs;
    
    if (problem_data->DomainVec().size() != 0)
    {
        // domain's material (2D or 3D)
        cmesh_v->SetDefaultOrder(problem_data->VelpOrder());
        
        cmesh_v->SetDimModel(dimension);
        
        cmesh_v->SetAllCreateFunctionsContinuous();
        
        auto *mat = new TPZNullMaterial<>(problem_data->DomainVec()[0].matID);
        mat->SetNStateVariables(dimension);
        cmesh_v->InsertMaterialObject(mat);
        
        materialIDs.insert(problem_data->DomainVec()[0].matID);
        
        // boundary condition's material
        TPZFMatrix<STATE> val1(1, 1, 0.0);
        TPZManVector<STATE> val2(1, 0.0);
        
        for (const auto &bc : problem_data->TangentialBCs())
        {
            val2 = bc.value;
            auto BCmat = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }
        
        cmesh_v->AutoBuild(materialIDs);
        
        int64_t nCon = cmesh_v->NConnects();
        for (int64_t i = 0; i < nCon; i++)
        {
            TPZConnect &newnod = cmesh_v->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        gmesh->ResetReference();
    }

    cmesh_v->ExpandSolution();
    
    problem_data->MeshVector()[0] = cmesh_v;
    
    return cmesh_v;
}
// *        END          *

// **********************
//    PRESS COMP MESH
// **********************
TPZCompMesh *CreateCMeshP(ProblemData *problem_data, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("CMesh_p");
    
    const int dimension = problem_data->Dim();
    std::set<int> materialIDs;
    
    if (problem_data->DomainVec().size() != 0)
    {
        
        cmesh_p->SetDimModel(dimension);
        
        cmesh_p->SetDefaultOrder(problem_data->VelpOrder() - 1);
        cmesh_p->SetAllCreateFunctionsContinuous();
        
        // domain's material
        auto *mat = new TPZNullMaterial<>(problem_data->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);
        
        materialIDs.insert(problem_data->DomainVec()[0].matID);
        
        cmesh_p->AutoBuild(materialIDs);
        
        int64_t nCon = cmesh_p->NConnects();
        for (int64_t i = 0; i < nCon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(2);
        }
        
        gmesh->ResetReference();
        
        materialIDs.clear();
    }
    
    cmesh_p->ExpandSolution();
    
    problem_data->MeshVector()[1] = cmesh_p;
    
    return cmesh_p;
}
// *        END          *


// **********************
// MULTPHYSICS COMP MESH
// **********************
TPZMultiphysicsCompMesh *CreateMultiphysicsCMesh(ProblemData *problem_data, TPZGeoMesh *gmesh)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);
    
    cmesh_m->SetName("CMesh_M");
    
    cmesh_m->SetDefaultOrder(problem_data->VelpOrder());
    cmesh_m->SetAllCreateFunctionsContinuous();
    
    // creating materials
    if (problem_data->DomainVec().size() != 0)
    {
        const int dimension = problem_data->Dim();
        STATE viscosity = problem_data->DomainVec()[0].viscosity;
        
        // 1. for domain
        TPZStokesMaterialTH *mat = new TPZStokesMaterialTH(problem_data->DomainVec()[0].matID, dimension, viscosity);
        
        cmesh_m->InsertMaterialObject(mat);
        
        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.0);
        TPZManVector<STATE> val2(3, 0.0);
        
        for (const auto &bc : problem_data->TangentialBCs())
        {
            val2 = bc.value;
            
            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            
            if (mat->HasExactSol())
                matBC2->SetForcingFunctionBC(mat->ExactSol(), 4);
            
            cmesh_m->InsertMaterialObject(matBC);
        }
    }
    
    TPZManVector<int, 2> active_approx_spaces(problem_data->MeshVector().size(), 1);
    
    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, problem_data->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();
    
    return cmesh_m;
}
// *        END          *


// **********************
//    SOLVER FUNCTION
// **********************
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m)
{
    TPZFStructMatrix<STATE> matskl(cmesh_m);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    // setting direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU // ECholesky // ELDLt
    an.SetSolver(step);
    
    // assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    // solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}
// *        END          *


// **********************
//    PRINT RESULTS
// **********************
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name)
{
    std::cout << "--------- PostProcess ---------" << std::endl;
    TPZSimpleTimer postProc("Post Process time");
    
    TPZVec<std::string> fields;
    
    if (problem_data->HasAnalyticSolution())
        fields = {
            "Pressure",
            "ExactPressure",
            "ErrorPressure",
            
            "Velocity",
            "ExactVelocity",
            "ErrorVelocity",
            
            "Stress",
            "ExactStress",
            "ErrorStress"
        };
    
    else
     fields = {
            "Pressure",
        
            "Velocity",
        
            "Stress",
        };
    
    auto vtk = TPZVTKGenerator(cmesh_m, fields, file_name, problem_data->Resolution());
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}
// *        END          *


// **********************
//        ERROS
// **********************
void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh, ProblemData *problem_data)
{
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
    an.SetThreadsForError(global_nthread);
    
    TPZMaterial *mat = cmesh->FindMaterial(problem_data->DomainVec()[0].matID);
    TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(mat);
    TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
    
    bool store_errors = false;
    std::ofstream ErrorOut(file_name + "_Errors.txt");
    
    TPZSimpleTimer time_error;
    
    an.PostProcessError(Errors, store_errors);
    std::cout << "Error time = " << time_error.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    ErrorOut << "###### Computed Errors ######" << std::endl;
    for (int i = 0; i < Errors.size(); i++)
        ErrorOut << "L2 " << fields[i] << " = " << Errors[i] << std::endl;
}
// *        END          *
