#include <iostream>
#include <cmath>
#include<fstream>
#include <string>

#include <pzgmesh.h>
#include <TPZGenGrid2D.h>
#include <TPZVTKGeoMesh.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <pzfmatrix.h>
#include <TPZGeoElement.h>
#include <TPZInterfaceEl.h>
#include <TPZMultiphysicsInterfaceEl.h>
#include <TPZGeoLinear.h>
#include <TPZNullMaterial.h>

#include "TPZMeshOperator.h"
#include "ProblemData.h"

TPZGeoMesh* CreateGMesh(SimulationData* simData){
    TPZGeoMesh* gmesh = new TPZGeoMesh;

    int ndiv = simData->nDiv();
    TPZVec<int> divisions = {ndiv, ndiv, ndiv};
    
    TPZGenGrid2D grid(divisions , simData->X0(), simData->X1());
    grid.SetElementType(MMeshType::EQuadrilateral);

    grid.Read(gmesh);
    grid.SetBC(gmesh, 4, -1);
    grid.SetBC(gmesh, 5, -4);
    grid.SetBC(gmesh, 6, -2);
    grid.SetBC(gmesh, 7, -3);

    InsertInterfaceElement(gmesh);

    gmesh -> BuildConnectivity();

    return gmesh;
}

void InsertInterfaceElement(TPZGeoMesh* gmesh){
    int64_t nEl = gmesh -> NElements();
    
    for(int64_t el = 0; el < nEl; el++){
        TPZGeoEl* geoEl = gmesh->Element(el);
        
        if (!geoEl) continue;
        if (geoEl->HasSubElement()) continue;
        
        int nside = geoEl->NSides();
        
        for(int side = 0; side < nside; side++){
            if(geoEl->SideDimension(side) != gmesh->Dimension()-1) continue;
            
            TPZGeoElSide geoElSide(geoEl, side);
            TPZGeoElSide neighbour = geoElSide.Neighbour();
            
            if (neighbour == geoElSide) continue;
            if(neighbour.Element()->HasSubElement()) continue;
            
            while (neighbour != geoElSide){
                
                if(neighbour.Element()->Dimension() == gmesh->Dimension()-1){
                    int neighbourMatId = neighbour.Element()->MaterialId();
                    
                    if(neighbourMatId == -1){
                        TPZGeoElBC(geoElSide, 11);
                        
                    } else if (neighbourMatId == -2){
                        TPZGeoElBC(geoElSide, 12);
                        
                    } else if(neighbourMatId == -3){
                        TPZGeoElBC(geoElSide, 13);
                        
                    } else if(neighbourMatId == -4){
                        TPZGeoElBC(geoElSide, 14);
                        
                    }
                    break;
                }
    
                neighbour = neighbour.Neighbour();
            }
            
            if(neighbour == geoElSide) TPZGeoElBC(geoElSide, 4);
        }
    }
}

TPZCompMesh* CreateCMeshV(SimulationData* simData, TPZGeoMesh* gmesh){
    TPZCompMesh* cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v -> SetDefaultOrder(simData -> VelpOrder());
    cmesh_v->SetDimModel(simData->DomainVec()[0].dim); // change i in the case of more than 1 domain

    // domain's material - 2D
    auto* mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_v -> InsertMaterialObject(mat);

    cmesh_v -> SetAllCreateFunctionsHDiv();

    // boundary conditions' material
    TPZFMatrix<STATE> val1(1.,1.,0.);
    TPZManVector<STATE> val2(2.,0.);

    auto BCbott = mat -> CreateBC(mat, -1, 0, val1, val2);
    cmesh_v -> InsertMaterialObject(BCbott);

    auto BCtop = mat -> CreateBC(mat, -2, 0, val1, val2);
    cmesh_v -> InsertMaterialObject(BCtop);

    auto BCleft = mat -> CreateBC(mat, -3, 0, val1, val2);
    cmesh_v -> InsertMaterialObject(BCleft);

    auto BCright = mat -> CreateBC(mat, -4, 0, val1, val2);
    cmesh_v -> InsertMaterialObject(BCright);
    
    int64_t ncel = cmesh_v -> NElements();
    
    for(int64_t i = 0; i < ncel; i++){
        TPZCompEl* compEl = cmesh_v -> ElementVec()[i];
        if(!compEl) continue;
        
        TPZInterfaceElement* interfaceEl = dynamic_cast<TPZInterfaceElement*>(compEl);
        if(interfaceEl) DebugStop();
    }

    cmesh_v -> AutoBuild();
    cmesh_v -> AdjustBoundaryElements();
    cmesh_v -> CleanUpUnconnectedNodes();
    
    return cmesh_v;
}

TPZCompMesh* CreateCmeshP(SimulationData* simData, TPZGeoMesh* gmesh){
    TPZCompMesh* cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p -> SetDefaultOrder(simData->VelpOrder());
    cmesh_p -> SetDimModel(simData->DomainVec()[0].dim);
    
    cmesh_p -> SetAllCreateFunctionsContinuous();
    cmesh_p -> ApproxSpace().CreateDisconnectedElements(true);
    
    // domain's material
    auto* mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_p -> InsertMaterialObject(mat);
    
    // traction on interface element material
    auto matLambda = new TPZNullMaterial<>(4);
    cmesh_p -> InsertMaterialObject(matLambda);
    
    // traction on boundary material
    auto matLambdaBC_bott = new TPZNullMaterial<>(11);
    cmesh_p -> InsertMaterialObject(matLambdaBC_bott);
    
    auto matLambdaBC_top = new TPZNullMaterial<>(12);
    cmesh_p -> InsertMaterialObject(matLambdaBC_top);
    
    auto matLambdaBC_left = new TPZNullMaterial<>(13);
    cmesh_p -> InsertMaterialObject(matLambdaBC_left);
    
    auto matLambdaBC_right = new TPZNullMaterial<>(14);
    cmesh_p -> InsertMaterialObject(matLambdaBC_right);
    
    int64_t ncel = cmesh_p->NElements();
    for(int64_t i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh_p->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
    }
    
    std::set<int> materialIDs;
    materialIDs.insert(simData -> DomainVec()[0].matID);
    
    cmesh_p -> AutoBuild(materialIDs);
    gmesh -> ResetReference();
    
    materialIDs.clear();
    materialIDs.insert(4);
    materialIDs.insert(11);
    materialIDs.insert(12);
    materialIDs.insert(13);
    materialIDs.insert(14);
    
    cmesh_p -> SetAllCreateFunctionsContinuous();
    cmesh_p -> ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh_p->SetDefaultOrder(simData->TracpOrder());
    cmesh_p->SetDimModel(simData->DomainVec()[0].dim - 1);
    cmesh_p->AutoBuild(materialIDs);
    
    int64_t ncon = cmesh_p -> NConnects();
    for(int64_t i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh_p -> ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh_p -> ExpandSolution();
    
    return cmesh_p;
}

void PrintMesh(TPZGeoMesh* gmesh, TPZCompMesh* cmesh_v, TPZCompMesh* cmesh_p){
    std::ofstream VTKGeoMeshFile("GMeshtest.vtk");
    std::ofstream TextGeoMeshFile("GMeshtest.txt");
    
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, VTKGeoMeshFile);
    gmesh->Print(TextGeoMeshFile);
    
    std::ofstream VTKCMeshVFile("CMesh_V.vtk");
    std::ofstream TextCMeshVFile("CMesh_V.txt");
    
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_v, VTKCMeshVFile);
    cmesh_v->Print(TextCMeshVFile);
    
    std::ofstream VTKCMeshPFile("CMesh_P.vtk");
    std::ofstream TextCMeshPFile("CMesh_P.txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_p, VTKCMeshPFile);
    cmesh_p->Print(TextCMeshPFile);
}