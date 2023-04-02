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
#include <TPZNullMaterialCS.h>

#include "TPZStokesMaterial.h"
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
    TPZManVector<STATE> val2(1, 0.);
    
    for(const auto& bc : simData->BCs()){
        if(bc.type == 0 || bc.type == 1){
            val2[0] = bc.value;
            
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

TPZMultiphysicsCompMesh* TPZMeshOperator::CreateMultiPhysicsMesh(ProblemData* simData, TPZGeoMesh* gmesh){
    TPZMultiphysicsCompMesh* cmesh_m = new TPZMultiphysicsCompMesh(gmesh);
    cmesh_m->SetName("MultiPhysics Mesh - Stokes Material");
    cmesh_m->SetDefaultOrder(simData->VelpOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    //Creating Materials
    //1. For domain
    TPZStokesMaterial* material = new TPZStokesMaterial(simData->DomainVec()[0].matID,simData->Dim());
    cmesh_m->InsertMaterialObject(material);

    // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);

    for(const auto& bc : simData->BCs()){
        val2[0] = bc.value;

        TPZBndCond* matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
        cmesh_m->InsertMaterialObject(matBC);
    }
    
    // 2.1 - Material for 1D tangential traction
//    TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->InterfaceID());
//    matLambda->SetDimension(simData->Dim()-1);
//    matLambda->SetNStateVariables(1);
//    cmesh_m->InsertMaterialObject(matLambda);

    // 2.2 - Material for interfaces (Inner)
//    TPZStokesMaterial *matInterfaceLeft = new TPZStokesMaterial(fmatInterfaceLeft,fdim);
//    matInterfaceLeft->SetSimulationData(f_sim_data);
//    matInterfaceLeft->SetMultiplier(1.);
//    cmesh->InsertMaterialObject(matInterfaceLeft);
//
//    TPZMHMNavierStokesMaterial *matInterfaceRight = new TPZMHMNavierStokesMaterial(fmatInterfaceRight,fdim);
//    matInterfaceRight->SetSimulationData(f_sim_data);
//    matInterfaceRight->SetMultiplier(-1.);
//    cmesh->InsertMaterialObject(matInterfaceRight);
//
//    // material para fixar uma pressao
//    TPZL2ProjectionCS<STATE> *l2proj = new TPZL2ProjectionCS<>(fmatPoint,fdim);
//    l2proj->SetNeedSol(1);
//    l2proj->SetScaleFactor(l2proj->BigNumber());
//    cmesh->InsertMaterialObject(l2proj);
    
    // insert an L2 projection element to the gmesh and the dist flux mesh
//    if(f_mesh_vector[3]->MaterialVec().find(fmatPoint) == f_mesh_vector[3]->MaterialVec().end())
//    {
//        // insert a volumetric element to gmesh
//        TPZCompMesh *pressmesh = f_mesh_vector[3];
//        int64_t nelc = pressmesh->NElements();
//        TPZCompEl *firstcel = 0;
//        TPZGeoEl *firstgel = 0;
//        for (int el = 0; el<nelc; el++) {
//            auto *cel = pressmesh->Element(el);
//            if(cel && cel->Reference() && cel->Reference()->Dimension() == fdim) {
//                firstcel = cel;
//                firstgel = cel->Reference();
//                break;
//            }
//        }
//        if(!firstgel) DebugStop();
//        TPZGeoElSide gelside(firstgel);
//        TPZGeoElBC gbc(gelside,fmatPoint);
//        TPZGeoEl *fixav = gbc.CreatedElement();
//        // insert an element in the pressure mesh
//        // first insert the material
//        TPZNullMaterial<> *fixavmat = new TPZNullMaterial<>(fmatPoint,fdim);
//        pressmesh->InsertMaterialObject(fixavmat);
//        // set the element type to discontinuous
//        pressmesh->SetAllCreateFunctionsDiscontinuous();
//        pressmesh->SetDimModel(fdim);
//        pressmesh->SetDefaultOrder(0);
//        auto discel = pressmesh->CreateCompEl(fixav);
//        pressmesh->ExpandSolution();
//        int64_t cindex = firstcel->ConnectIndex(0);
//        discel->Connect(0).DecrementElConnected();
//        discel->SetConnectIndex(0, cindex);
//        // adjust the solution vector
//        pressmesh->CleanUpUnconnectedNodes();
//    }
    
    // 3.1 - Material for 1D tangential traction on the boundaries
    //    val2[0] = 0.0;
//    TPZBndCond *matLambdaBC_bott = material->CreateBC(material, fmatLambdaBC_bott, fdirichlet_sigma, val1, val2);
//    cmesh->InsertMaterialObject(matLambdaBC_bott);
//
//    val2[0] = 0.;
//    TPZBndCond *matLambdaBC_top = material->CreateBC(material, fmatLambdaBC_top, fneumann_sigma, val1, val2);
//    cmesh->InsertMaterialObject(matLambdaBC_top);
//
//    val2[0] = 0.0;
//    TPZBndCond *matLambdaBC_left = material->CreateBC(material, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
//    cmesh->InsertMaterialObject(matLambdaBC_left);
//
//    TPZBndCond *matLambdaBC_right = material->CreateBC(material, fmatLambdaBC_right, fdirichlet_sigma, val1, val2);
//    cmesh->InsertMaterialObject(matLambdaBC_right);

//    int ncel = cmesh_m->NElements();
//    for(int i =0; i<ncel; i++){
//        TPZCompEl * compEl = cmesh->ElementVec()[i];
//        if(!compEl) continue;
//        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//        if(facel)DebugStop();
//
//    }
//
//    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
//
//    TPZManVector<int,5> active_approx_spaces(f_mesh_vector.size(),1);
//
//    cmesh->BuildMultiphysicsSpace(active_approx_spaces,f_mesh_vector);
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    return cmesh;
//
}

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
