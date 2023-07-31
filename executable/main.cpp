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

int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif
  
    bool printdata = false;

    std::string filepath = "../examples/";
    std::string filename = "AxisymmetricObstructedAxialFlow";

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
        TPZMeshOperator::CondenseElements(cmesh_m);
    }

    TPZLinearAnalysis an(cmesh_m, true);
    //TPZSSpStructMatrix<> strmat(cmesh_m);
    TPZFStructMatrix<> strmat(cmesh_m);
    //TPZSkylineStructMatrix<> strmat(cmesh_m);

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

    {
        TPZSimpleTimer timer("Solving", true);
        
        an.Assemble();
        
        an.Solve();
    }

    an.Solve();

    

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
    TPZVTKGenerator vtk(cmesh_m, {"Pressure", "Velocity", "Tension"}, filename, simData.Resolution());
    vtk.Do();

    std::cout << "\n\nSimulation finished without errors :) \n\n";
            
	return 0;
}
