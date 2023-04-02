#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"
#include "TPZStokesMaterial.h"

int main(){
    bool printmesh = false;

    std::string filenamejson =  "StokesData.json";

    ProblemData simData;
    simData.ReadJson(filenamejson);
    simData.Print();
    
    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    TPZManVector<TPZCompMesh*, 2> meshVector = {cmesh_v, cmesh_p};
    simData.SetMeshVector(meshVector);
    
//    TPZMultiphysicsCompMesh* cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);
    
    if(printmesh){
        TPZMeshOperator::PrintMesh(gmesh, cmesh_v, cmesh_p);
    }

    std::cout << "\n\nSimulation finished without errors :) \n\n";
	return 0;
}
