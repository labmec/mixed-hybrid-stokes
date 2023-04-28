#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"


int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("carlos.cfg");
#endif
    
    bool printmesh = true;

    std::string filenamejson =  "StokesData.json";

    ProblemData simData;
    simData.ReadJson(filenamejson);
//    simData.Print();

    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);

    TPZMultiphysicsCompMesh* cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);

    if(printmesh){
        TPZMeshOperator::PrintMesh(gmesh, cmesh_v, cmesh_p, cmesh_m);
    }

    TPZLinearAnalysis an(cmesh_m, true);
    TPZFStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(0);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    an.Assemble();
    
    an.Solve();

    std::cout << "\n\nSimulation finished without errors :) \n\n";
	return 0;
}
