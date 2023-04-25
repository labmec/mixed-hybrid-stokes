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
#include <pzintel.h>

#include "TPZInterfaceMaterial.h"
#include "TPZStokesMaterial.h"
#include "TPZMeshOperator.h"
#include "ProblemData.h"

TPZGeoMesh* TPZMeshOperator::CreateGMesh(ProblemData* simData){
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gmesh->SetName("Geometric Mesh");
   
    TPZGmshReader reader;
    
    reader.GeometricGmshMesh(simData->MeshName(), gmesh);

   TPZMeshOperator::InsertLambdaGEl(simData, gmesh);

    gmesh->BuildConnectivity();

    return gmesh;
}

void TPZMeshOperator::InsertLambdaGEl(ProblemData* simData, TPZGeoMesh* gmesh){
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
            
            if(neighbour == geoElSide) TPZGeoElBC(geoElSide, simData->LambdaID());
        }
    }
}

void TPZMeshOperator::InsertBCInterfaces(TPZMultiphysicsCompMesh* cmesh_m, ProblemData* simData, TPZGeoMesh* gmesh){
    TPZManVector<int, 2> Interfaces(2,0);
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();
    
    if(!gmesh) DebugStop();
    
    int64_t nel = gmesh->NElements();
    
    TPZVec<int> IDVec(simData->TracBCs().size(), 0);
    for(int i = 0; i<simData->TracBCs().size(); i++){
        IDVec[i] = simData->TracBCs()[i].matID;
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
                    // Check if it is working in the case with refined meshes
                    DebugStop();
                    
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
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();
    
    int dim = cmesh_m->Dimension();
    
    if(!gmesh) DebugStop();
    
    int nInterfaceCreated = 0;
    
    int matfrom = simData->LambdaID();
    
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
//    cmesh_v->SetDefaultOrder(simData->TracpOrder());
    cmesh_v->SetDefaultOrder(2);
    cmesh_v->SetDimModel(simData->Dim());

    cmesh_v->SetAllCreateFunctionsHDiv();
    
    // domain's material - 2D
    auto* mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_v->InsertMaterialObject(mat);
    
    // boundary conditions' material
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0.);
    
    for(const auto& bc : simData->VelBCs()){
            val2 = bc.value;
            
            auto BCmat = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
    }
    
    cmesh_v->AutoBuild();
    cmesh_v->LoadReferences();
    
    // setting the approximation order for the volume elements
    int64_t ncEl = cmesh_v->NElements();
    for(int64_t cEl=0; cEl<ncEl; cEl++){
        TPZCompEl* compEl = cmesh_v->Element(cEl);
        
        // only in those elements whose dimension equals to the simulation dim
        if(compEl->Dimension()==simData->Dim()){
            // dynamica casting the compEl object to use the ForceSideOrder function
            TPZInterpolatedElement* intercEl = dynamic_cast<TPZInterpolatedElement*>(compEl);
            
            // checking if the dynamic cast exists
            if(!intercEl) continue;
            
            // finally using the desired function
//            intercEl->ForceSideOrder(compEl->Reference()->NSides()-1, simData->VelpOrder()); 
        }
    }
    
    cmesh_v->CleanUpUnconnectedNodes();
    
    // expanding the solution vector
    cmesh_v->ExpandSolution();
    
    simData->MeshVector().Resize(2);
    simData->MeshVector()[0] = cmesh_v;
    
    return cmesh_v;
}

TPZCompMesh* TPZMeshOperator::CreateCmeshP(ProblemData* simData, TPZGeoMesh* gmesh){
    TPZCompMesh* cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("H1 - Pressure");
    cmesh_p->SetDefaultOrder(simData->VelpOrder());
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
    auto matLambda = new TPZNullMaterial<>(simData->LambdaID());
    cmesh_p->InsertMaterialObject(matLambda);
    
    materialIDs.insert(simData->LambdaID());
    
    // traction on boundary material
    for(const auto& bc : simData->TracBCs()){
            auto matLambdaBC = new TPZNullMaterial<>(bc.matID);
            cmesh_p->InsertMaterialObject(matLambdaBC);
            
            materialIDs.insert(bc.matID);
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
    TPZStokesMaterial* material = new TPZStokesMaterial(simData->DomainVec()[0].matID,simData->Dim(), simData->DomainVec()[0].viscosity);
    cmesh_m->InsertMaterialObject(material);

    // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,1.);

    for(const auto& bc : simData->VelBCs()){
        val2 = bc.value;

        TPZBndCond* matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
        cmesh_m->InsertMaterialObject(matBC);
    }
    
    for(const auto& bc : simData->TracBCs()){
        val2 = bc.value;

        TPZBndCond* matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
        cmesh_m->InsertMaterialObject(matBC);
    }
    
//     2.1 - Material for 1D tangential traction
    TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->LambdaID());
    matLambda->SetDimension(simData->Dim()-1);
    matLambda->SetNStateVariables(1);
    cmesh_m->InsertMaterialObject(matLambda);

    // 2.2 - Material for interfaces (Inner)
    TPZInterfaceMaterial *matInterfaceLeft = new TPZInterfaceMaterial(simData->InterfaceID(),simData->Dim());
    matInterfaceLeft->SetMultiplier(1.);
    cmesh_m->InsertMaterialObject(matInterfaceLeft);

    TPZInterfaceMaterial *matInterfaceRight = new TPZInterfaceMaterial(-simData->InterfaceID(), simData->Dim());
    matInterfaceRight->SetMultiplier(-1.);
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
