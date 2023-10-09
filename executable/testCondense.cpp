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
#include <pzmatred.h>

int main()
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif
  
    bool printdata = false;

    std::string filepath = "../examples/";
    std::string filename = "AxisymmetricHagenPoiseuilleFlow";

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

    cmesh_m->SaddlePermute();
    TPZLinearAnalysis an(cmesh_m,RenumType::ENone);

    TPZMeshOperator::PrintCompMesh(cmesh_m);

    int64_t nEq = cmesh_m->NEquations();
    int64_t nEq1d = 0;
    for (const TPZGeoEl* gel : gmesh->ElementVec())
    {
        int matid = gel->MaterialId();
        if (matid == simData.InterfaceID() || matid == simData.LambdaID() || matid == simData.DomainVec()[0].matID) //we only want the 1d elements
            continue;
        const TPZCompEl* cel = gel->Reference();
        if (cel)
        {
            const int nconnects = cel->NConnects();
            for (int i = 0; i < nconnects; i++)
            {
                TPZConnect con = cel->Connect(i);
                
            }
        }

    }
    TPZMatRed<STATE, TPZFMatrix<STATE>> *matRed = new TPZMatRed<STATE, TPZFMatrix<STATE>>(nEq,);

    TPZSSpStructMatrix<> strmat(cmesh_m);
    //TPZFStructMatrix<> strmat(cmesh_m);
    //TPZSkylineStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(0);

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
