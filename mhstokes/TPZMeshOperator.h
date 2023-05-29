#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"

class TPZMeshOperator{
public:
    static TPZGeoMesh* CreateGMesh(ProblemData* simData);
    static void InsertLambdaGEl(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCMeshV(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMv(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMp(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZMultiphysicsCompMesh* CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void CondenseElements(TPZMultiphysicsCompMesh* cmesh_m);
    
    static void InsertBCInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void InsertInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void PrintGeoMesh(TPZGeoMesh* gmesh);

    static void PrintCompMesh(TPZVec<TPZCompMesh*> CMeshVec);
};
#endif
