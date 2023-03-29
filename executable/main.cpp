#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"

int main(){
    bool printmesh = true;

//    std::string filenamejson =  "StokesData.json";
//
//    SimulationData simData;
//    simData.ReadJson(filenamejson);
//
//    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);
//    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);
//    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
//    
//    TPZMultiphysicsCompMesh* cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);
//
//    if(printmesh){
//        TPZMeshOperator::PrintMesh(gmesh, cmesh_v, cmesh_p);
//    }


    bool hdiv = false;

    if(hdiv){
        Example::HdivConforming();
    } else {
        Example::H1Conforming();
    }
    
    std::cout << "Simulation finishes without errors :) \n";
	return 0;
}
