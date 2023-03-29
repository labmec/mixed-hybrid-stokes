#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"

class TPZMeshOperator{
public:
    static TPZGeoMesh* CreateGMesh(SimulationData* simData);
    static void InsertInterfaceElement(TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCMeshV(SimulationData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshP(SimulationData* simData, TPZGeoMesh* gmesh);
    static void PrintMesh(TPZGeoMesh* gmesh, TPZCompMesh* cmesh_v, TPZCompMesh* cmesh_p);
    
    static TPZMultiphysicsCompMesh* CreateMultiPhysicsMesh(SimulationData* simData, TPZGeoMesh* gmesh);
};
#endif
