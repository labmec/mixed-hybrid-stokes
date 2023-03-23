#include <iostream>
#include <fstream>
#include <string>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"

int main(){
    bool printmesh = true;

    std::string filenamejson =  "StokesData.json";

    SimulationData simData;
    simData.ReadJson(filenamejson);
    
    TPZGeoMesh* gmesh = CreateGMesh(&simData);
    TPZCompMesh* cmesh_v = CreateCMeshV(&simData, gmesh);
    TPZCompMesh* cmesh_p = CreateCmeshP(&simData, gmesh);

    if(printmesh){
        PrintMesh(gmesh, cmesh_v, cmesh_p);
    }

	return 0;
}
