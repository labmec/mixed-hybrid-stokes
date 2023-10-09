#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"
#include <TPZAnalyticSolution.h>

class TPZMeshOperator{
public:
    enum SimSpace {EStandard, EConstant};
    
    static TPZGeoMesh* CreateGMesh(ProblemData* simData);
    
    static void InsertLambdaGEl(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCMeshV(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMv(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshMp(ProblemData* simData, TPZGeoMesh* gmesh);
        
    static TPZMultiphysicsCompMesh* CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh* cmesh_m);
    
    static void InsertBCInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void InsertInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void PrintGeoMesh(TPZGeoMesh* gmesh);

    static void PrintCompMesh(TPZVec<TPZCompMesh*> CMeshVec);
    
    static void PrintCompMesh(TPZCompMesh* cmesh);
    
    static void GenerateMshFile(ProblemData* geoFile);

    static void CheckSideOrientOfCompEl(ProblemData* simData, TPZGeoMesh* gmesh);
};
#endif
