#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"
#include <TPZAnalyticSolution.h>

class TPZMeshOperator{
public:
    
    static TPZGeoMesh* CreateGMesh(ProblemData* simData);
    
    static void InsertLagrangeMultipliers(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCMeshV(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshG(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static TPZCompMesh* CreateCmeshPm(ProblemData* simData, TPZGeoMesh* gmesh);
        
    static TPZMultiphysicsCompMesh* CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh, TPZAnalyticSolution* sol);
    
    static void CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh* cmesh_m, TPZGeoMesh *gmesh);
    
    static void InsertInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void PrintGeoMesh(TPZGeoMesh* gmesh);

    static void PrintCompMesh(TPZVec<TPZCompMesh*> CMeshVec);
    
    static void PrintCompMesh(TPZCompMesh* cmesh);
    
    static void GenerateMshFile(ProblemData* geoFile);

    static void CheckSideOrientOfCompEl(ProblemData* simData, TPZGeoMesh* gmesh);
    
    static void ConfigureObstructionFilter(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, std::set<int64_t> &removeEquations);
    
    static void ConfigureBoundaryFilter(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, std::set<int64_t> &removeEquations);
};
#endif
