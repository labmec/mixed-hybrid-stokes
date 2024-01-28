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
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif
  
    bool printdata = false;

    std::string filepath = "../DataInput/";
    std::string filename = "DifferentObstructions";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = simData.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::ENone;
    
    // creating the meshes
    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);

    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    if(simData.CondensedElements() && simData.HdivType() != ProblemData::EConstant){
        TPZCompMesh* cmesh_Mp = TPZMeshOperator::CreateCmeshPm(&simData, gmesh);
        TPZCompMesh* cmesh_Mv = TPZMeshOperator::CreateCmeshG(&simData, gmesh);
    }

    TPZMultiphysicsCompMesh *cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh, flow);
    
    if (simData.CondensedElements())
        TPZMeshOperator::CondenseElements(&simData, cmesh_m, gmesh);
    
    cmesh_m->SaddlePermute();

    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    TPZLinearAnalysis an(cmesh_m, RenumType::ENone);
    
//    TPZSSpStructMatrix<STATE, TPZStructMatrixOT<STATE>> strmat(cmesh_m);
    TPZSSpStructMatrix<> strmat(cmesh_m);
//    TPZFStructMatrix<> strmat(cmesh_m);
//     TPZSkylineStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(0);
    
    // Applying filters to null normal and tangential velocities
    TPZEquationFilter filter(cmesh_m->NEquations());
    std::set<int64_t> removeEquations;
    
    // on the boundary
    TPZMeshOperator::ConfigureBoundaryFilter(gmesh, cmesh_m, &simData, removeEquations);
    
    // between the domain and the obstruction
    if (simData.ObstructionID() != - 1)
        TPZMeshOperator::ConfigureObstructionFilter(gmesh, cmesh_m, &simData, removeEquations);
    
    filter.ExcludeEquations(removeEquations);
    strmat.EquationFilter() = filter;

    TPZSimpleTimer timer;
    
    an.SetStructuralMatrix(strmat);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    std::cout << "Starting assemble...\n";
    an.Assemble();
    std::cout << "Finished assemble...\n";
    
    std::cout << "Starting solver...\n";
    std::cout << "Number of Equations " << cmesh_m->NEquations() << std::endl;

    TPZFMatrix<REAL> aux_solution(cmesh_m->NEquations(), 1);
    aux_solution.Zero();
    
    an.Solve();
    std::cout << "Finished solver...\n";
    
    std::cout << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << std::endl;

    // vtk export
    TPZVTKGenerator vtk(cmesh_m, {"Pressure", "Velocity", "Stress"}, filename, simData.Resolution());
    vtk.Do();
    
    std::cout << "\n\nSimulation finished without errors :) \n\n";
       
	return 0;
}
