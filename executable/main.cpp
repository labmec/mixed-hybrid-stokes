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
  
    bool printdata = true;

    std::string filepath = "../examples/Elasticity/";
    std::string filename = "UniformTension3D";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");

    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);

    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);

    if(simData.CondensedElements()){
        TPZCompMesh* cmesh_Mp = TPZMeshOperator::CreateCmeshMp(&simData, gmesh);
        TPZCompMesh* cmesh_Mv = TPZMeshOperator::CreateCmeshMv(&simData, gmesh);
    }

    TPZMultiphysicsCompMesh *cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);

    if (simData.CondensedElements())
    {
        TPZMeshOperator::CondenseElements(&simData, cmesh_m);
    }
    TPZMeshOperator::PrintCompMesh(cmesh_m);
    //cmesh_m->SaddlePermute();

    TPZLinearAnalysis an(cmesh_m, RenumType::ENone);
    
    //TPZSSpStructMatrix<STATE, TPZStructMatrixOT<STATE>> strmat(cmesh_m);
    TPZSSpStructMatrix<> strmat(cmesh_m);
    //TPZFStructMatrix<> strmat(cmesh_m);
    // TPZSkylineStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(0);

    // TPZEquationFilter filter(cmesh_m->NEquations());
    // std::set<int64_t> setremove;

    // for (auto el : cmesh_m->ElementVec())
    // {
    //     auto matid = el->Reference()->MaterialId();
    //     if (matid != 6)
    //         continue;
    //     int64_t nconnects = el->NConnects();
    //     for (int64_t ic = 0; ic < nconnects; ic++)
    //     {
    //         TPZConnect &c = el->Connect(ic);
    //         int64_t blocknumber = c.SequenceNumber();
    //         auto firsteq = cmesh_m->Block().Position(blocknumber);
    //         int64_t blocksize = cmesh_m->Block().Size(blocknumber);
    //         for (int64_t eq = firsteq; eq < firsteq + blocksize; eq++)
    //         {
    //             setremove.insert(eq);
    //         }
    //     }
    // }
    // filter.ExcludeEquations(setremove);
    // strmat.EquationFilter() = filter;
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    std::cout << "Starting assemble...\n";
    an.Assemble();
    std::cout << "Finished assemble...\n";

    std::cout << "Starting solver...\n";
    an.Solve();
    std::cout << "Finished solver...\n";

    if (printdata)
    {
        simData.Print();

        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        TPZFMatrix<STATE> &Sol=an.Solution();
        Sol.Print("Sol=", std::cout , EMathematicaInput);
    }

    TPZManVector<REAL,10> Errors(3);
    //an.PostProcessError(Errors, false);

    // vtk export
    TPZVTKGenerator vtk(cmesh_m, {"Pressure", "Displacement", "Stress", "Strain"}, filename, simData.Resolution());
    vtk.SetNThreads(0);
    vtk.Do();

    std::cout << "\n\nSimulation finished without errors :) \n\n";
       
	return 0;
}
