#ifndef STOKESEXAMPLE_H
#define STOKESEXAMPLE_H

#include "ProblemData.h"

TPZGeoMesh* CreateGMesh(SimulationData* simData);
void InsertInterfaceElement(TPZGeoMesh* gmesh);

TPZCompMesh* CreateCMeshV(SimulationData* simData, TPZGeoMesh* gmesh);

TPZCompMesh* CreateCmeshP(SimulationData* simData, TPZGeoMesh* gmesh);
void PrintMesh(TPZGeoMesh* gmesh, TPZCompMesh* cmesh_v, TPZCompMesh* cmesh_p);

#endif
