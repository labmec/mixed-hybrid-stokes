#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include <fstream>
#include <sstream>
#include <string>

#include "TPZAnalyticSolution.h"
#include "TPZInterfaceAxisymStokesMaterial.h"
#include "TPZInterface1dStokesMaterial.h"
#include "TPZInterface1dFlux.h"
#include "TPZStokesMaterial.h"
#include "TPZAxisymStokesMaterial.h"
#include "TPZ1dStokesMaterial.h"
#include "TPZMixedLinearElasticMaterial.h"
#include "ProblemData.h"

class TPZMeshOperator{
public:
    
    static TPZGeoMesh* CreateGMesh(ProblemData* simData, bool blend);
    
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
    
//    static void SetExactArcRepresentation(const TPZAutoPointer<TPZGeoMesh> &gmesh, ProblemData *simData);
    static void SetExactArcRepresentation(TPZGeoMesh &gmesh, ProblemData *simData);
    
//    static void SetExactCylinderRepresentation(const TPZAutoPointer<TPZGeoMesh> &gmesh, ProblemData *simData);
    static void SetExactCylinderRepresentation(TPZGeoMesh &gmesh, ProblemData *simData);
    
    static void printVTKWJacInfo(std::string filename, TPZGeoMesh *gmesh);
};
#endif
