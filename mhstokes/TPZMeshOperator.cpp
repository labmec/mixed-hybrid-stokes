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
#include <TPZMultiphysicsCompMesh.h>
#include <TPZGeoLinear.h>
#include <TPZNullMaterial.h>
#include <TPZGMshReader.h>

#include "TPZNavierStokesMaterial.h"
#include "TPZMeshOperator.h"
#include "ProblemData.h"

TPZGeoMesh* TPZMeshOperator::CreateGMesh(ProblemData* simData){
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gmesh->SetName("Geometric Mesh");
   
    TPZGmshReader reader;
    
    reader.GeometricGmshMesh(simData->MeshName(), gmesh);

   TPZMeshOperator::InsertInterfaceElement(simData, gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

void TPZMeshOperator::InsertInterfaceElement(ProblemData* simData, TPZGeoMesh* gmesh){
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
                    
                    break;
                }
    
                neighbour = neighbour.Neighbour();
            }
            
            if(neighbour == geoElSide) TPZGeoElBC(geoElSide, simData->InterfaceID());
        }
    }
}

TPZCompMesh* TPZMeshOperator::CreateCMeshV(ProblemData* simData, TPZGeoMesh* gmesh){
    TPZCompMesh* cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v->SetName("Hdiv Mesh - Velocity");
    cmesh_v->SetDefaultOrder(simData->VelpOrder());
    cmesh_v->SetDimModel(simData->Dim());

    cmesh_v -> SetAllCreateFunctionsHDiv();
    
    // domain's material - 2D
    auto* mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_v -> InsertMaterialObject(mat);

    // boundary conditions' material
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 1.);
    
    for(const auto& bc : simData->BCs()){
        if(bc.type == 0 || bc.type == 1){
            val2[0] *= bc.value;
            
            auto BCmat = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
        }
    }
    
    cmesh_v -> AutoBuild();
    cmesh_v -> CleanUpUnconnectedNodes();
    
    return cmesh_v;
}

TPZCompMesh* TPZMeshOperator::CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh){
    TPZCompMesh* cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("H1 - Pressure");
    cmesh_p->SetDefaultOrder(simData->TracpOrder());
    cmesh_p->SetDimModel(simData->Dim());
    
    cmesh_p->SetAllCreateFunctionsContinuous();
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
    
    // domain's material
    auto* mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_p->InsertMaterialObject(mat);
    
    std::set<int> materialIDs;
    materialIDs.insert(simData->DomainVec()[0].matID);
    
    cmesh_p->AutoBuild(materialIDs);
    gmesh->ResetReference();
    
    materialIDs.clear();
    
    // matlambda traction material
    auto matLambda = new TPZNullMaterial<>(simData->InterfaceID());
    cmesh_p->InsertMaterialObject(matLambda);
    
    materialIDs.insert(simData->InterfaceID());
    
    // traction on boundary material
    for(const auto& bc : simData->BCs()){
        if(bc.type == 4 || bc.type == 3){
            auto matLambdaBC = new TPZNullMaterial<>(bc.matID);
            cmesh_p->InsertMaterialObject(matLambdaBC);
            
            materialIDs.insert(bc.matID);
        }
    }
    
    cmesh_p->SetAllCreateFunctionsContinuous();
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh_p->SetDefaultOrder(simData->TracpOrder());
    cmesh_p->SetDimModel(simData->Dim() - 1);
    cmesh_p->AutoBuild(materialIDs);
    
    int64_t ncon = cmesh_p->NConnects();
    for(int64_t i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh_p->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh_p->ExpandSolution();
    
    return cmesh_p;
}

//TPZMultiphysicsCompMesh* TPZMeshOperator::CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh){
//    TPZMultiphysicsCompMesh* cmesh_m = new TPZMultiphysicsCompMesh(gmesh);
//    cmesh_m->SetName("MultiPhysics Mesh - Stokes Material");
//    cmesh_m->SetDefaultOrder(simData->VelpOrder());
//    cmesh_m->SetAllCreateFunctionsMultiphysicElem();
//
//    //Creating Materials
//    //1. For domain
//    TPZNavierStokesMaterial* material = new TPZNavierStokesMaterial(simData->DomainVec()[0].matID,simData->Dim());
//    cmesh_m->InsertMaterialObject(material);
//    
//    
//    
//
//}

void TPZMeshOperator::PrintMesh(TPZGeoMesh* gmesh, TPZCompMesh* cmesh_v, TPZCompMesh* cmesh_p){
    std::cout << "\nPrinting geometric and computational meshes in .txt and .vtk formats...";
    
    std::ofstream VTKGeoMeshFile("GMesh.vtk");
    std::ofstream TextGeoMeshFile("GMesh.txt");
    
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
