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

   TPZMeshOperator::InsertInterfaceMaterial(simData, gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

void TPZMeshOperator::InsertInterfaceMaterial(ProblemData* simData, TPZGeoMesh* gmesh){
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

void TPZMeshOperator::InsertBCInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh){
    TPZManVector<int, 2> Interfaces(2,0);
    Interfaces[0] = simData->LambdaID();
    Interfaces[1] = -simData->LambdaID();
    
    if(!gmesh) DebugStop();
    
    int64_t nel = gmesh->NElements();
    
    TPZVec<int> IDVec(simData->BCs().size(), 0);
    for(int i = 0; i<simData->BCs().size(); i++){
        IDVec[i] = simData->BCs()[i].matID;
    }
    
    for(auto const& BcMatID : IDVec){
        for(int64_t el=0; el<nel; el++){
            TPZGeoEl* gel = gmesh->Element(el);
            int meshDim = gmesh->Dimension();
            int matID = gel->MaterialId();
            
            if(matID!=BcMatID) continue;
            
            int nsides = gel->NSides();
            TPZGeoElSide gelSide(gel, nsides-1);
            TPZCompElSide celSide = gelSide.Reference();
            
            TPZStack<TPZGeoElSide> neighbourSet;
            gelSide.AllNeighbours(neighbourSet);
            
            int64_t nneighs = neighbourSet.size();
            
            TPZManVector<int64_t, 1> LeftElIndex(1,0), RightElIndex(1,0);
            LeftElIndex[0] = 0;
            RightElIndex[0] = 1;
            
            for(int stack_i = 0; stack_i < nneighs; stack_i++){
                TPZGeoElSide neigh = neighbourSet[stack_i];
                int neighMatID = neigh.Element()->MaterialId();
                TPZCompElSide celNeigh = neigh.Reference();
                
                int64_t neighIndex = neigh.Element()->Index();
                
                if(neigh.Element()->Dimension()!=meshDim) continue;
                
                if(neigh.Element()->HasSubElement()){
                    DebugStop();
                    // Check if it is working in the case with refined meshes
//                    TPZStack<TPZGeoElSide> subElSide;
//                    neigh.GetAllSiblings(subElSide);
//
//                    for(int i_sub = 0; i_sub < subElSide.size(); i_sub++){
//                        TPZCompElSide celSubNeigh = subElSide[i_sub].Reference();
//                        TPZGeoElBC gbcSub(subElSide[i_sub], BcMatID);
//
//                        TPZMultiphysicsInterfaceElement* interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbcSub.CreatedElement(), celSubNeigh, celSide);
//                        interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
//                    }
                    
                } else {
                    
                    TPZGeoElBC gbc(gelSide, Interfaces[0]);
                    
                    TPZMultiphysicsInterfaceElement* interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                    interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                }
            }
        }
    }
}

void TPZMeshOperator::InsertInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh){
    TPZManVector<int, 2> Interfaces(2,0);
    Interfaces[0] = simData->LambdaID();
    Interfaces[1] = -simData->LambdaID();
    
    int dim = cmesh_m->Dimension();
    
    if(!gmesh) DebugStop();
    
    int nInterfaceCreated = 0;
    
    int matfrom = simData->InterfaceID();
    
    int64_t nel = gmesh->NElements();
    
    for(int64_t el=0; el<nel; el++){
        TPZGeoEl* gel = gmesh->Element(el);
        int meshdim = gmesh->Dimension();
        int matid = gel->MaterialId();
        
        if(matid != matfrom) continue;
        if(gel->HasSubElement()) continue;
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel, nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        TPZStack<TPZGeoElSide> neighbourSet;
        gelside.AllNeighbours(neighbourSet);
        
        gelside.LowerLevelCompElementList2(1);
        
        int64_t nneighs = neighbourSet.size();
        
        TPZManVector<int64_t, 3> LeftElIndices(1,0.), RightElIndices(1,0);
        LeftElIndices[0] = 0;
        RightElIndices[0] = 1;
        
        for(int stack_i = 0; stack_i < nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourSet[stack_i];
            int neighMatID = neigh.Element()->MaterialId();
            TPZCompElSide celneigh = neigh.Reference();
            
            if(!celside) DebugStop();
            
            if(neigh.Element()->HasSubElement()){
                DebugStop();
//                TPZStack<TPZGeoElSide> subElements;
//
//                TPZStack<TPZGeoElSide> subEl;
//                neigh.GetAllSiblings(subEl);
//
//                for(int iSub = 0; iSub<subEl.size(); iSub++){
//                    TPZCompElSide celSubNeigh = subEl[iSub].Reference();
//                    TPZGeoElBC gbc_sub(subEl[iSub], Interfaces[stack_i]);
//
//                    TPZMultiphysicsInterfaceElement* interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc_sub.CreatedElement(), celSubNeigh, celside);
//                    interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
//
//                    nInterfaceCreated++;
//                }
                
            } else {
                
                int64_t neighIndex = neigh.Element()->Index();
                
                if(neigh.Element()->Dimension()!=meshdim) continue;
                
                TPZGeoElBC gbc(gelside, Interfaces[stack_i]);
                
                TPZMultiphysicsInterfaceElement* interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celneigh, celside);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;
            }
        }
    }
    
    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created " << nInterfaceCreated << std::endl;
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
    
    simData->MeshVector().Resize(2);
    simData->MeshVector()[0] = cmesh_v;
    
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
    
    simData->MeshVector()[1] = cmesh_p;
    
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
    
//     2.1 - Material for 1D tangential traction
    TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->InterfaceID());
    matLambda->SetDimension(simData->Dim()-1);
    matLambda->SetNStateVariables(1);
    cmesh_m->InsertMaterialObject(matLambda);

    // 2.2 - Material for interfaces (Inner)
    TPZStokesMaterial *matInterfaceLeft = new TPZStokesMaterial(simData->LambdaID(),simData->Dim());
    cmesh_m->InsertMaterialObject(matInterfaceLeft);

    TPZStokesMaterial *matInterfaceRight = new TPZStokesMaterial(-simData->LambdaID(),simData->Dim());
    cmesh_m->InsertMaterialObject(matInterfaceRight);

    // Creating computational elements that will manage the mesh approximation space:
    int64_t ncel = cmesh_m->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh_m->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();

    }

    TPZManVector<int,3> active_approx_spaces(simData->MeshVector().size(),1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces,simData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    TPZMeshOperator::InsertBCInterfaces(cmesh_m, simData, gmesh);
    TPZMeshOperator::InsertInterfaces(cmesh_m, simData, gmesh);

    return cmesh_m;

}

void TPZMeshOperator::PrintMesh(TPZGeoMesh* gmesh, TPZCompMesh* cmesh_v, TPZCompMesh* cmesh_p, TPZMultiphysicsCompMesh* cmesh_m){
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
    
    std::ofstream VTKCMeshMFile("CMesh_M.vtk");
    std::ofstream TextCMeshMFile("CMesh_M.txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, VTKCMeshMFile);
    cmesh_m->Print(TextCMeshMFile);
}
