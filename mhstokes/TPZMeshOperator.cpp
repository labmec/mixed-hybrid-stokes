#include <iostream>
#include <cmath>
#include <fstream>
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
#include <TPZNullMaterialCS.h>
#include <pzintel.h>
#include <pzelementgroup.h>
#include <pzcondensedcompel.h>
#include <filesystem>
#include <tpzchangeel.h>

#include "TPZMeshOperator.h"
#include "ProblemData.h"
#include "pzreal.h"

void TPZMeshOperator::GenerateMshFile(ProblemData *simData)
{
    std::string command;

    command = "gmsh " + simData->MeshName() + ".geo -o " + simData->MeshName() + ".msh -algo del2d -2 -order 1";
    system(command.c_str());
}

TPZGeoMesh *TPZMeshOperator::CreateGMesh(ProblemData *simData)
{

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetName("GeoMesh");

    if (simData->CreateMshFile())
    {
        GenerateMshFile(simData);
    }

    TPZGmshReader reader;

    reader.GeometricGmshMesh(simData->MeshName() + ".msh", gmesh);
    
    //     using geo blend
    if (simData->CsvFile() != "none")
    {
        switch (simData->Dim()) {
            case 2:
                {
                    simData->ReadCirclesData();
                    TPZMeshOperator::SetExactArcRepresentation(*gmesh, simData);
                }
                break;
                
            case 3:
            {
                simData->ReadCylindersData();
                TPZMeshOperator::SetExactCylinderRepresentation(*gmesh, simData);
            }
            default:
                break;
        }
    }
    
    TPZMeshOperator::InsertLagrangeMultipliers(simData, gmesh);
    
    gmesh->BuildConnectivity();

    return gmesh;
}

void TPZMeshOperator::InsertLagrangeMultipliers(ProblemData *simData, TPZGeoMesh *gmesh)
{
    int64_t nEl = gmesh->NElements();
    
    // we look for tangential BCs
    TPZVec<int> IDVec(simData->TangentialBCs().size(), 0);
    
    for (int i=0; i<simData->TangentialBCs().size(); i++)
        IDVec[i] = simData->TangentialBCs()[i].matID;
    
    for (auto const &BcMatID : IDVec)
    {
        for (int64_t el = 0; el < nEl; el++)
        {
            int meshDim = gmesh->Dimension();
            
            TPZGeoEl *geoEl = gmesh->Element(el);
            int matID = geoEl->MaterialId();
            
            if (matID != BcMatID)
                continue;
            
            int nSides = geoEl->NSides();
            TPZGeoElSide geoElSide(geoEl, nSides-1);
            TPZCompElSide compElSide = geoElSide.Reference();
            
            TPZStack<TPZGeoElSide> neighbourSet;
            geoElSide.AllNeighbours(neighbourSet);
            
            int64_t nneighs = neighbourSet.size();
            
            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neighbour = neighbourSet[stack_i];
                int neighMatID = neighbour.Element()->MaterialId();
                TPZCompElSide compElNeigh = neighbour.Reference();
                
                int64_t neighIndex = neighbour.Element()->Index();
                
                if (neighbour.Element()->Dimension() != meshDim)
                    continue;
                
                if(neighbour.Element()->HasSubElement())
                    DebugStop();
                
                TPZGeoElBC(neighbour, simData->InterfaceID());
                
                neighbour = neighbour.Neighbour();
                
                if (neighbour.Element()->MaterialId() != simData->InterfaceID())
                    DebugStop();
                
                TPZGeoElBC(neighbour, simData->LambdaID());
                
                neighbour = neighbour.Neighbour();
                
                if (neighbour.Element()->MaterialId() != simData->LambdaID())
                    DebugStop();
                
                TPZGeoElBC(neighbour, 21);
            }
        }
    }
    
    // we look for two domain neighbour elements
    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);

        if (!geoEl)
            continue;
        if (geoEl->HasSubElement())
            continue;
        if (geoEl->Dimension() != gmesh->Dimension())
            continue;

        int nside = geoEl->NSides();

        if (simData->DomainVec().size() != 0)
        {
            for (int side = 0; side < nside; side++)
            {
                if (geoEl->SideDimension(side) != gmesh->Dimension() - 1)
                    continue;
                
                TPZGeoElSide geoElSide(geoEl, side);
                TPZGeoElSide neighbour = geoElSide.Neighbour();

                if (neighbour == geoElSide)
                    continue;
                if (neighbour.Element()->HasSubElement())
                    continue;

                TPZGeoElSide neighbour2 = neighbour;
                while (neighbour2 != geoElSide)
                {
                    if (neighbour2.Element()->MaterialId() == simData->InterfaceID())
                    {
                        break;
                    }
                    
                    if (neighbour2.Element()->Dimension() == gmesh->Dimension())
                        neighbour = neighbour2;
                        
                    neighbour2 = neighbour2.Neighbour();
                }

                if (neighbour2 == geoElSide)
                {
                    TPZGeoElBC(geoElSide, simData->InterfaceID());
                    TPZGeoElBC(neighbour, simData->InterfaceID());
                    
                    neighbour2 = neighbour.Neighbour();
                    neighbour = geoElSide.Neighbour();
                    
                    if(neighbour.Element()->MaterialId() != simData->InterfaceID() || neighbour2.Element()->MaterialId() != simData->InterfaceID())
                        DebugStop();
                    
                    TPZGeoElBC(neighbour, simData->LambdaID());
                    TPZGeoElBC(neighbour2, simData->LambdaID());
                    
                    neighbour = neighbour.Neighbour();
                    neighbour2 = neighbour2.Neighbour();
                    
                    if (neighbour.Element()->MaterialId() != simData->LambdaID() || neighbour2.Element()->MaterialId() != simData->LambdaID())
                        DebugStop();
                    
                    TPZGeoElBC(neighbour, 21);
                    TPZGeoElBC(neighbour2, 21);
                    
                    neighbour = neighbour.Neighbour();
                    
                    if (neighbour.Element()->MaterialId() != 21)
                        DebugStop();
                    
                    TPZGeoElBC(neighbour, simData->TanVelID());
                }
            }
        }
    }
    
    gmesh->BuildConnectivity();
}

void TPZMeshOperator::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh)
{
    if (!gmesh)
        DebugStop();
    
    TPZManVector<int64_t, 3> LeftElIndices(1,0), RightElIndices(1,1);
    
    int nInterfaceCreated = 0;
    
    int mat_lambda = simData->LambdaID();
    int mat_tan = simData->TanVelID();
    
    int meshDim = cmesh_m->Dimension();
    
    std::set<int> BcIDs;
    for (int i = 0; i < simData->TangentialBCs().size(); i++)
        BcIDs.insert(simData->TangentialBCs()[i].matID);
    
    int64_t nEl = gmesh->NElements();
    
    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);
        int matID = geoEl->MaterialId();
        
        if (geoEl->Dimension() != meshDim)
            continue;
        
        if (geoEl->HasSubElement())
            continue;
        
        int nSides = geoEl->NSides();
        
        for (int side = 0; side < nSides; side++)
        {
            if (geoEl->SideDimension(side) != meshDim-1)
                continue;
            
            TPZGeoElSide geoElSide(geoEl, side);
            TPZGeoElSide geoEl_neighbour = geoElSide.Neighbour();
            int neighbour_matID = geoEl_neighbour.Element()->MaterialId();
            
            if (neighbour_matID != simData->InterfaceID())
            {
                int a = 0;
            }
            
            if (neighbour_matID == simData->InterfaceID())
            {
                TPZCompElSide compElSide = geoElSide.Reference();
                
                if (!compElSide)
                    DebugStop();
                
                TPZGeoElSide geoEl_neighbour2 = geoEl_neighbour.Neighbour();
                
                if (geoEl_neighbour2.Element()->MaterialId() == simData->ObstructionID())
                    DebugStop();
                
                if (geoEl_neighbour2.Element()->MaterialId() != simData->LambdaID())
                    continue;
                
                TPZCompElSide compEl_neighbour = geoEl_neighbour2.Reference();
                
                if (!compEl_neighbour)
                    DebugStop();
                
                if (geoEl_neighbour2.Element()->HasSubElement())
                    DebugStop();
                
                // creating the interface between the Hdiv domain element and tangential stress
                TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, geoEl_neighbour.Element(), compElSide, compEl_neighbour);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;
                
                // we take the second interface between tangential stress and tangential velocity
                geoEl_neighbour = geoEl_neighbour2.Neighbour();
                
                if (geoEl_neighbour.Element()->MaterialId() != 21)
                    DebugStop();
                
                TPZGeoElSide geoEl_neighbour3 = geoEl_neighbour;
                for (; geoEl_neighbour3 != geoElSide; geoEl_neighbour3++)
                {
                    int neighbour_matID2 = geoEl_neighbour3.Element()->MaterialId();
                                      
                    if (neighbour_matID2 != mat_tan && BcIDs.find(neighbour_matID2) == BcIDs.end())
                        continue;
                    
                    if (neighbour_matID2 == simData->ObstructionID())
                        DebugStop();
                    
                    compElSide = geoEl_neighbour3.Reference();
                    
                    if (!compElSide)
                        DebugStop();
                    
                    if (geoEl_neighbour3.Element()->HasSubElement())
                        DebugStop();
                    
                    
                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, geoEl_neighbour.Element(), compElSide, compEl_neighbour);
                    interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                    nInterfaceCreated++;
                    
                    break;
                }
            }
        }
    }
    
    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created: " << nInterfaceCreated << std::endl;
}

TPZCompMesh *TPZMeshOperator::CreateCMeshV(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v->SetName("Cmesh_V");
    
    std::set<int> materialID;
    
    cmesh_v->ApproxSpace().CreateDisconnectedElements(false);
    
    // domain's material (2D or 3D)
    if (simData->DomainVec().size() != 0)
    {
        cmesh_v->SetDefaultOrder(simData->VelpOrder());
        cmesh_v->SetDimModel(simData->Dim());
        
        if (simData->HdivType() == ProblemData::EConstant)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        }
        else if (simData->HdivType() == ProblemData::EStandard)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
        }
        
        cmesh_v->SetAllCreateFunctionsHDiv();
        
        auto *mat_normal = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_v->InsertMaterialObject(mat_normal);
        
        materialID.insert(simData->DomainVec()[0].matID);
        
        // boundary condition's material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);
        
        for (const auto &bc: simData->NormalBCs())
        {
            val2 = bc.value;
            
            auto BCmat = mat_normal->CreateBC(mat_normal, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
            materialID.insert(bc.matID);
        }
        
        cmesh_v->AutoBuild(materialID);
        
        // Increasing internal function order
        int64_t ncEl = cmesh_v->NElements();
        
        for (int64_t cEl = 0; cEl < ncEl; cEl++)
        {
            TPZCompEl *compEl = cmesh_v->Element(cEl);
            
            if (compEl->Dimension() == simData->Dim())
            {
                TPZInterpolatedElement *intercEl = dynamic_cast<TPZInterpolatedElement *>(compEl);
                
                // checking of the dynamic cast existis
                if (!intercEl)
                    continue;
                
                intercEl->ForceSideOrder(compEl->Reference()->NSides()-1, simData->VelpOrder()+2);
            }
        }
        
        int64_t ncon = cmesh_v->NConnects();
        for (int64_t i =0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_v->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(2);
        }
        
        gmesh->ResetReference();
        materialID.clear();
        
        // tangential velocity material
        auto mat_tan = new TPZNullMaterial<>(simData->TanVelID());
        mat_tan->SetNStateVariables(simData->Dim()-1); // in 3d, there are 2 state variables (one at each tangential direction)
        cmesh_v->InsertMaterialObject(mat_tan);
        
        materialID.insert(simData->TanVelID());
        
        // tangential velocity on boundary material
        for (const auto &bc : simData->TangentialBCs())
        {
            auto matBC = new TPZNullMaterial<>(bc.matID);
            matBC->SetNStateVariables(simData->Dim()-1);
            cmesh_v->InsertMaterialObject(matBC);
            
            materialID.insert(bc.matID);
        }
        
        cmesh_v->ApproxSpace().CreateDisconnectedElements(true);
        cmesh_v->SetDefaultOrder(simData->TracpOrder());
        cmesh_v->SetDimModel(simData->Dim()-1);
        cmesh_v->AutoBuild(materialID);
        
        ncon = cmesh_v->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_v->ConnectVec()[i];
            if (newnod.LagrangeMultiplier() == 0)
                newnod.SetLagrangeMultiplier(3);
        }
        
        gmesh->ResetReference();
    }
    
    cmesh_v->ExpandSolution();
    simData->MeshVector()[0] = cmesh_v;
    
    return cmesh_v;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshP(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("CMesh_P");
    
    std::set<int> materialID;
    
    if (simData->DomainVec().size() != 0)
    {
        cmesh_p->SetDimModel(simData->Dim());
        
        if (simData->HdivType() == ProblemData::EConstant)
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->SetDefaultOrder(0);
        }
        else if (simData->HdivType() == ProblemData::EStandard)
        {
            cmesh_p->SetDefaultOrder(simData->VelpOrder()+2);
            cmesh_p->SetAllCreateFunctionsContinuous();
        }
        
        cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        
        // domain's material
        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);
        
        materialID.insert(simData->DomainVec()[0].matID);
        
        cmesh_p->AutoBuild(materialID);
        gmesh->ResetReference();
        
        materialID.clear();
        
        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(4);
        }
        
        // matLambda traction material
        auto matLambda = new TPZNullMaterial<>(simData->LambdaID());
        matLambda->SetNStateVariables(simData->Dim()-1);
        cmesh_p->InsertMaterialObject(matLambda);
        
        materialID.insert(simData->LambdaID());
        
        if (simData->TracpOrder()>0)
        {
            cmesh_p->SetAllCreateFunctionsContinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        
        cmesh_p->SetDefaultOrder(simData->TracpOrder());
        cmesh_p->SetDimModel(simData->Dim()-1);
        cmesh_p->AutoBuild(materialID);
        
        gmesh->ResetReference();
        
        ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            
            if (newnod.LagrangeMultiplier() == 0)
                newnod.SetLagrangeMultiplier(1);
        }
    }
    
    cmesh_p->ExpandSolution();
    simData->MeshVector()[1] = cmesh_p;
    
    return cmesh_p;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshG(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_g = new TPZCompMesh(gmesh);
    
    cmesh_g->SetName("CMesh_g");
    cmesh_g->SetDefaultOrder(0);
    cmesh_g->SetDimModel(simData->Dim());
    
    cmesh_g->SetAllCreateFunctionsDiscontinuous();
    
    auto mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_g->InsertMaterialObject(mat);
    
    cmesh_g->AutoBuild();
    cmesh_g->AdjustBoundaryElements();
    cmesh_g->CleanUpUnconnectedNodes();
    
    int64_t ncon = cmesh_g->NConnects();
    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &newnod = cmesh_g->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(5);
    }
    
    simData->MeshVector()[2] = cmesh_g;
    
    return cmesh_g;
}

TPZCompMesh *TPZMeshOperator::CreateCmeshPm(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_pm = new TPZCompMesh(gmesh);
    
    cmesh_pm->SetName("CMesh_pm");
    cmesh_pm->SetDefaultOrder(0);
    cmesh_pm->SetDimModel(simData->Dim());
    
    cmesh_pm->SetAllCreateFunctionsDiscontinuous();
    
    auto mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_pm->InsertMaterialObject(mat);
    
    cmesh_pm->AutoBuild();
    cmesh_pm->AdjustBoundaryElements();
    cmesh_pm->CleanUpUnconnectedNodes();
    
    int64_t ncon = cmesh_pm->NElements();
    
    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &newnod = cmesh_pm->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(5);
    }
    
    simData->MeshVector()[3] = cmesh_pm;
    
    return cmesh_pm;
}

TPZMultiphysicsCompMesh *TPZMeshOperator::CreateMultiPhysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh, TPZAnalyticSolution *sol)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");

    cmesh_m->SetDefaultOrder(simData->VelpOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials
    std::set<int> materialID;
    
    // 1. For domain
    if (simData->DomainVec().size() != 0)
    {
        REAL viscosity = simData->DomainVec()[0].viscosity;
        
        if (dynamic_cast<TStokesAnalytic*>(sol))
        {
            TStokesAnalytic *flow = dynamic_cast<TStokesAnalytic*>(sol);
            viscosity = flow->fvisco;
            
            if (flow->fExactSol == TStokesAnalytic::ENone)
                sol = nullptr;
        }
        
        TPZStokesMaterial *material = new TPZStokesMaterial(simData->DomainVec()[0].matID, simData->Dim(), viscosity);
        
        if (sol) material->SetExactSol(sol->ExactSolution(), 3);
        if (sol) material->SetForcingFunction(sol->ForceFunc(), 3);
            
        cmesh_m->InsertMaterialObject(material);
        materialID.insert(simData->DomainVec()[0].matID);
        
        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3,3,0.);
        TPZManVector<STATE> val2(3,0.);
        
        // Normal
        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;
            
            TPZBndCond *matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
            
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            if (sol) matBC2->SetForcingFunctionBC(sol->ExactSolution(), 3);
            
            cmesh_m->InsertMaterialObject(matBC);
            materialID.insert(bc.matID);
        }
        
        // Tangential
        for (const auto &bc : simData->TangentialBCs())
        {
            val2 = bc.value;
            
            TPZBndCond *matBC = material->CreateBC(material, bc.matID, bc.type, val1, val2);
            
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            if (sol) matBC2->SetForcingFunctionBC(sol->ExactSolution(), 3);
            
            cmesh_m->InsertMaterialObject(matBC);
            materialID.insert(bc.matID);
        }
        
        // 3. Material for Tangential Velocity
        TPZNullMaterialCS<> *mat_tan = new TPZNullMaterialCS<>(simData->TanVelID());
        mat_tan->SetDimension(simData->Dim()-1);
        mat_tan->SetNStateVariables(simData->Dim()-1);
        
        cmesh_m->InsertMaterialObject(mat_tan);
        materialID.insert(simData->TanVelID());
        
        // 4. Material for Tangential Traction
        TPZNullMaterialCS<> *mat_lambda = new TPZNullMaterialCS<>(simData->LambdaID());
        mat_lambda->SetDimension(simData->Dim()-1);
        mat_lambda->SetNStateVariables(simData->Dim()-1);
        
        cmesh_m->InsertMaterialObject(mat_lambda);
        materialID.insert(simData->LambdaID());
    }
    
    TPZManVector<int, 2> active_approx_spaces(simData->MeshVector().size(), 1);
    
    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, simData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();
    
    // 5. Material for Interfaces (inner)
    TPZInterfaceAxisymStokesMaterial *matInterface = new TPZInterfaceAxisymStokesMaterial(simData->InterfaceID(), simData->Dim()-1);
    matInterface->SetMultiplier(1.);
    
    cmesh_m->InsertMaterialObject(matInterface);
    
    TPZInterfaceAxisymStokesMaterial *matInterface2 = new TPZInterfaceAxisymStokesMaterial(21, simData->Dim()-1);
    matInterface2->SetMultiplier(1.);
    
    cmesh_m->InsertMaterialObject(matInterface2);
    
    InsertInterfaces(cmesh_m, simData, gmesh);
    
    if (simData->CondensedElements())
    {
        cmesh_m->SetName("CMesh_M_BeforeCond");
        cmesh_m->ComputeNodElCon();
    }
    
    return cmesh_m;
}

void TPZMeshOperator::CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh *cmesh_m, TPZGeoMesh *gemesh)
{
    bool condensedPressure = (simData->MeshVector().size()==2 && simData->HdivType() != ProblemData::EConstant) ? false : true;
    int64_t nCompEl = cmesh_m->ElementVec().NElements();
    int dim = gemesh->Dimension();
    
    std::set<int64_t> externalNodes;
    std::vector<int64_t> groupIndex;
    groupIndex.reserve(nCompEl);
    
    TPZStack<TPZElementGroup *> elGroups;
    int count = 0;
    
    std::set<int> BCsIDs;
    for (int i = 0; i < simData->TangentialBCs().size(); i++)
        BCsIDs.insert(simData->TangentialBCs()[i].matID);
    for (int i = 0; i < simData->NormalBCs().size(); i++)
        BCsIDs.insert(simData->NormalBCs()[i].matID);
    
    // Creating the element groups for the domain
    for (int64_t el = 0; el < nCompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);
        
        if (compEl->Dimension() != dim)
            continue;
        
        int nConnect = compEl->NConnects();
        
        if (condensedPressure)
        {
            int64_t coIndex = compEl->ConnectIndex(nConnect - 1);
            externalNodes.insert(coIndex);
        }
        
        count++;
        groupIndex.push_back(compEl->Index());
        
        TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
        elGroups.Push(groupEl);
        elGroups[count-1]->AddElement(compEl);
    }
    
    int64_t nGeoEl = gemesh->NElements();
    for (int elgr = 0; elgr < elGroups.size(); elgr++)
    {
        TPZElementGroup *groupEl = elGroups[elgr];
        TPZCompEl *compEl = groupEl->GetElGroup()[0];
        
        if (!compEl) 
            DebugStop();
        
        TPZGeoEl *geoEl = compEl->Reference();
        
        if (geoEl->Dimension() != dim)
            DebugStop();
        
        if (geoEl->HasSubElement())
            DebugStop();
        
        int nSides = geoEl->NSides();
        
        int64_t compEl_Index = compEl->Index();
        
        for (int side = 0; side < nSides; side++)
        {
            if (geoEl->SideDimension(side) != dim - 1)
                continue;
            
            TPZGeoElSide geoEl_Side(geoEl, side);
            
            TPZStack<TPZGeoElSide> allNeighbours;
            geoEl_Side.AllNeighbours(allNeighbours);
            
            TPZStack<TPZCompEl*> neighbourCompEls;
            
            if (allNeighbours.size() > 0)
            {
                
                if (allNeighbours[0].Element()->MaterialId() == simData->InterfaceID())
                {
                    for (int i = 0; i < 3; i++)
                    {
                        TPZGeoEl *gel = allNeighbours[i].Element();
                        int matID = gel->MaterialId();
                        
#ifdef PZDEBUG
                        if (matID != simData->InterfaceID() && matID != simData->LambdaID() && matID != 21)
                            DebugStop();
#endif
                        TPZCompEl *cel = allNeighbours[i].Element()->Reference();
                        
                        if (!cel) 
                            DebugStop();
                        
                        groupEl->AddElement(cel);
                    }
                }
            }
        }
    }
    
    cmesh_m->ComputeNodElCon();
    
    for (auto it = externalNodes.begin(); it != externalNodes.end(); it++)
    {
        int64_t coIndex = *it;
        cmesh_m->ConnectVec()[coIndex].IncrementElConnected();
    }
    
    // Creating condensed elements
    int64_t nEnvel = elGroups.NElements();
    for (int64_t iEnv = 0; iEnv < nEnvel; iEnv++)
    {
        TPZElementGroup *elGroup = elGroups[iEnv];
        new TPZCondensedCompElT<STATE>(elGroup);
    }
    
    cmesh_m->SetName("CMesh_M_Condensed");
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
}

void TPZMeshOperator::PrintGeoMesh(TPZGeoMesh *gmesh)
{
    std::cout << "\nPrinting geometric mesh in .txt and .vtk formats...\n";

    std::ofstream VTKGeoMeshFile(gmesh->Name() + ".vtk");
    std::ofstream TextGeoMeshFile(gmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, VTKGeoMeshFile);
    gmesh->Print(TextGeoMeshFile);
}

void TPZMeshOperator::PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting multiphysics mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}

void TPZMeshOperator::PrintCompMesh(TPZVec<TPZCompMesh *> CMeshVec)
{
    std::cout << "\nPrinting computational meshes in .txt and .vtk formats...\n";

    for (auto &cmesh : CMeshVec)
    {

        std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
        std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
        cmesh->ComputeNodElCon();
        cmesh->Print(TextCompMeshFile);
    }
}

void TPZMeshOperator::CheckSideOrientOfCompEl(ProblemData* simData, TPZGeoMesh* gmesh)
{
    if (simData->AxisymmetryDomainVec().size() != 0)
    {
        std::set<int> bcIds;
        for(auto& bc : simData->AxisymmetryBCs())
            bcIds.insert(bc.matID);

        int matid1d = simData->AxisymmetryDomainVec()[0].matID;

        for(TPZGeoEl* gel : gmesh->ElementVec())
        {
            const int gelmatid = gel->MaterialId();

            if (gelmatid != matid1d) //Check only 1d elements
                continue;

            int nside = gel->NSides(0); //Check only the vertices
            
            TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(gel->Reference());
            if (!intel) DebugStop();

            for (int is = 0; is < nside; is++)
            {
                const int sideorientgel = intel->GetSideOrient(is);
                TPZGeoElSide gelside(gel,is);
                //TPZGeoElSide neig = gelside.HasNeighbour(matid1d);
                TPZGeoElSide neig = gelside.Neighbour();

                for (; neig != gelside; neig++)
                {
                    int neigmatid = neig.Element()->MaterialId();

                    if (neigmatid == matid1d) //vertex is has an 1d element as neighbour
                    {
                        TPZInterpolatedElement* intelneig = dynamic_cast<TPZInterpolatedElement*>(neig.Element()->Reference());
                        const int sideorientneig = intelneig->GetSideOrient(neig.Side());

                        if (sideorientgel * sideorientneig != -1)
                        {
                            intelneig->SetSideOrient(neig.Side(), -sideorientgel);
                        }
                    }
                    else if (bcIds.find(neigmatid) != bcIds.end()) //vertex is has a bc as neighbour
                    {
                        intel->SetSideOrient(is, 1);
                    }
                }
            }
        }
    }
}

void TPZMeshOperator::ConfigureObstructionFilter(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, std::set<int64_t> &removeEquations)
{
    // Filter Function
    std::set<int64_t> removeConnectsSeq;
    
    for (auto el : gmesh->ElementVec())
    {
        auto matid = el->MaterialId();
        
        if (matid != simData->ObstructionID())
            continue;
        
        int64_t nSide = el->NSides();
        TPZGeoElSide geoElSide (el, nSide - 1);
        
        TPZStack<TPZGeoElSide> neighbourSet;
        geoElSide.AllNeighbours(neighbourSet);
        
        // Filtering tangential velocity
        for (auto neighbour : neighbourSet)
        {
            TPZGeoEl *geoEl_neighbour = neighbour.Element();
            
            if (geoEl_neighbour->MaterialId() != simData->TanVelID())
                continue;
            
            TPZCompEl *compEl_neighbour = geoEl_neighbour->Reference();
            
            int64_t nconnects = compEl_neighbour->NConnects();
            
            for (int64_t ic = 0; ic < nconnects; ic++)
                removeConnectsSeq.insert(compEl_neighbour->Connect(ic).SequenceNumber());
        }
        
        // Filtering normal velocity
        std::vector<TPZCompEl*> domain_neighbours;
        
        for (auto neighbour : neighbourSet)
        {
            TPZGeoEl *geo_neighbour = neighbour.Element();
            
            if (geo_neighbour->MaterialId() != simData->DomainVec()[0].matID)
                continue;
            
            domain_neighbours.push_back(geo_neighbour->Reference());
            
            if (domain_neighbours.size() == 2)
                break;
        }
        
        if (domain_neighbours.size() != 2)
            DebugStop();
        
        TPZCompEl *neighbour1 = domain_neighbours[0];
        TPZCompEl *neighbour2 = domain_neighbours[1];
        
        int64_t nConnects1 = neighbour1->NConnects();
        int64_t nConnects2 = neighbour2->NConnects();
        
        if (nConnects1 != nConnects2)
            DebugStop();
        
        std::set<int64_t> setIndex1;
        std::set<int64_t> setIndex2;
        
        neighbour1->BuildConnectList(setIndex1);
        neighbour2->BuildConnectList(setIndex2);
        
        std::set<int64_t> setIntersection;
        
        std::set_intersection(setIndex1.begin(), setIndex1.end(), setIndex2.begin(), setIndex2.end(), std::inserter(setIntersection, setIntersection.begin()));
        
        if (setIntersection.size() != 1)
            DebugStop();
        
        int64_t indexIntersec = *setIntersection.begin();
        
        removeConnectsSeq.insert(cmesh_m->ConnectVec()[indexIntersec].SequenceNumber());
    }
    
    for (auto blocknumber : removeConnectsSeq)
    {
        auto firsteq = cmesh_m->Block().Position(blocknumber);
        
        int64_t blocksize = cmesh_m->Block().Size(blocknumber);
        
        for (int64_t eq = firsteq; eq < firsteq + blocksize; eq++)
            removeEquations.insert(eq);
    }
}

void TPZMeshOperator::ConfigureBoundaryFilter(TPZGeoMesh *gmesh, TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, std::set<int64_t> &removeEquations)
{
    // Filter Function
    std::set<int64_t> removeConnectSeq;
    int NoSlipBC, NoPenetrationBC;
    
    for (auto bc : simData->TangentialBCs())
        if (bc.name == "NoSlip")
        {
            NoSlipBC = bc.matID;
            break;
        }
    
    for (auto bc : simData->NormalBCs())
        if (bc.name == "NoPenetration")
        {
            NoPenetrationBC = bc.matID;
            break;
        }
    
    for (auto el : gmesh->ElementVec())
    {        if (el->Dimension() != gmesh->Dimension() - 1)
            continue;
        
        auto elMatID = el->MaterialId();
        
        if (elMatID != NoSlipBC && elMatID != NoPenetrationBC)
            continue;

        TPZCompEl *compEl = el->Reference();
        
        int64_t nConnects = compEl->NConnects();
        
        if (nConnects != 1)
            DebugStop();
        
        removeConnectSeq.insert(compEl->Connect(0).SequenceNumber());
    }
    
    for (auto blockNumber : removeConnectSeq)
    {
        auto firstEq = cmesh_m->Block().Position(blockNumber);
        
        int64_t blockSize = cmesh_m->Block().Size(blockNumber);
        
        for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++)
            removeEquations.insert(eq);
    }
}

//void TPZMeshOperator::SetExactArcRepresentation(TPZAutoPointer<TPZGeoMesh> gmesh, ProblemData *simData)
void TPZMeshOperator::SetExactArcRepresentation(TPZGeoMesh &gmesh, ProblemData *simData)
{
    TPZVec<ProblemData::ArcData> circles = simData->ArcDataVec();

    std::map<int, int> arc_ids;
    std::map<int, bool> found_arcs;
    std::map<int, int> arc_newIDs;
    
    std::set<int> new_IDs;
    int availableID = 0;
    
    // calculate availableID
    for (auto domain : simData->DomainVec())
        availableID += domain.matID;

    for (auto bc : simData->NormalBCs())
        availableID += bc.matID;

    for (auto bc : simData->TangentialBCs())
        availableID += bc.matID;
    
    availableID *= 10;
    
    for (int i = 0; i < circles.size(); i++)
    {
        arc_ids[circles[i].matID] = i;
        arc_newIDs[circles[i].matID] = i + availableID;
        new_IDs.insert(i+availableID);
        found_arcs[circles[i].matID] = false;
    }
    
    TPZManVector<TPZGeoElSide, 50> neighs;
    
    for (auto el : gmesh.ElementVec())
    {
        // this way we avoid processing recently inserted elements
        if (!el || el->IsLinearMapping() == false) continue;
        
        const int matID = el->MaterialId();
        const bool is_arc = arc_ids.find(matID) != arc_ids.end();
        
        if (is_arc) // found arc
        {
            const int nSides = el->NSides();
            const int nNodes = el->NCornerNodes();
            
            found_arcs[matID] = true;
            const int arc_pos = arc_ids[matID];
            const REAL r = circles[arc_pos].radius;
            const REAL xc = circles[arc_pos].xc;
            const REAL yc = circles[arc_pos].yc;
            const REAL zc = circles[arc_pos].zc;
            
            auto new_el = el->CreateBCGeoEl(nSides - 1, arc_newIDs[matID]);
            auto *arc = TPZChangeEl::ChangeToArc3D(&gmesh, new_el->Index(), {xc, yc, zc}, r);
        
            // now we need to replace each neighbour by a blend element
            // so it can be deformed accordingly
            
            TPZGeoElSide geoElSide(arc, arc->NSides() - 1);
            
            // now we iterate through all the neighbours of the linear side
            TPZGeoElSide neighbour = geoElSide.Neighbour();
            
            // let us store all the neighbours
            std::set<TPZGeoElSide> all_neighbours;
            while (neighbour.Exists() && neighbour != geoElSide)
            {
                all_neighbours.insert(neighbour);
                neighbour = neighbour.Neighbour();
            }
            
            // let us replace all the neighbours
            for (auto neighbour : all_neighbours)
            {
                const auto neigh_side = neighbour.Side();
                auto neigh_el = neighbour.Element();
                
                /* 
                 let us take into account the possibility that
                 one triangle might be neighbour of two cylinders
                */
                
                if (!neigh_el->IsGeoBlendEl())
                {
                    const auto neigh_index = neigh_el->Index();
                    TPZChangeEl::ChangeToGeoBlend(&gmesh, neigh_index);
                }
                else
                {
                    neigh_el->SetNeighbourForBlending(neigh_side);
                }
            }
        }
    }
    
    for (auto arc : found_arcs)
    {
        if (!arc.second)
        {
            PZError<<__PRETTY_FUNCTION__
            <<"\n arc " << arc.first << " not found in mesh" << std::endl;
        }
    }
}

//void TPZMeshOperator::SetExactCylinderRepresentation(TPZAutoPointer<TPZGeoMesh> gmesh,  ProblemData *simData)
void TPZMeshOperator::SetExactCylinderRepresentation(TPZGeoMesh &gmesh,  ProblemData *simData)
{
    TPZVec<ProblemData::CylinderData> cylinders = simData->CylinderDataVec();

    std::map<int, int> cylinder_ids;
    std::map<int, bool> found_cylinders;
    std::map<int, int> cylinder_newIDs;

    int availableID = 0;
    std::set<int> new_IDs;

    // calculate availableID
    for (auto domain : simData->DomainVec())
        availableID += domain.matID;

    for (auto bc : simData->NormalBCs())
        availableID += bc.matID;

    for (auto bc : simData->TangentialBCs())
        availableID += bc.matID;

    availableID *= 10;
    
    for (auto i = 0; i < cylinders.size(); i++)
    {
        cylinder_ids[cylinders[i].matID] = i;
        cylinder_newIDs[cylinders[i].matID] = i + availableID;
        new_IDs.insert(i + availableID);
        found_cylinders[cylinders[i].matID] = false;
    }
    
    TPZManVector<TPZGeoElSide, 50> neighs;
    
    std::set<TPZGeoEl*> blend_neighs;
    for (auto el : gmesh.ElementVec())
    {
        // this way we avoid processing recently inserted elements
        if (!el || el->IsLinearMapping() == false) continue;
        
        const int matID = el->MaterialId();
        const bool is_cylinder = cylinder_ids.find(matID) != cylinder_ids.end();
        
        if (is_cylinder)
        {
            const int nSides = el->NSides();
            const int nNodes = el->NCornerNodes();
            
            found_cylinders[matID] = true;
            const int cylinder_pos = cylinder_ids[matID];
            const auto &cylData = cylinders[cylinder_pos];
            const REAL r = cylData.radius;
            TPZManVector<REAL, 3> xc = {cylData.xc, cylData.yc, cylData.zc};
            TPZManVector<REAL, 3> axis = {cylData.xaxis, cylData.yaxis, cylData.zaxis};
            
//            auto new_el = el->CreateBCGeoEl(nSides - 1, cylinder_newIDs[matID]);
            auto cyl = TPZChangeEl::ChangeToCylinder(&gmesh, el->Index(), xc, axis);
            
            // let us store all the neighbours
            std::set<TPZGeoElSide> all_neighbours;
            {
//                TPZGeoElSide elSide(el, nSides - 1);
//                all_neighbours.insert(elSide);
            }
            
            // just to ensure we wont remove any element twice
            std::set<TPZGeoEl*> neigh_els;
            
//            neigh_els.insert(el);
            
            for (auto is = cyl->NNodes(); is < cyl->NSides(); is++)
            {
                /*
                now we need to replace each neighbour by a blend element with a blend element 
                so it can be deformed accordingly
                */

                TPZGeoElSide geoElSide(cyl, is);

                // now we iterate through all the neighbours of the given side
                TPZGeoElSide neighbour = geoElSide.Neighbour();

                while (neighbour.Exists() && neighbour != geoElSide)
                {
                    auto neigh_el = neighbour.Element();

                    /*
                    let us skip:
                    1. neighbours that have been already identified (will be in neigh_els)
                    2. neighbours that are in the cylinder wall as well
                    */
        
                   auto check_el = neigh_els.find(neigh_el) == neigh_els.end();
                   if (check_el && neigh_el->MaterialId() != matID)
                   {
                       all_neighbours.insert(neighbour);
                       neigh_els.insert(neigh_el);
                   }
                   neighbour = neighbour.Neighbour();
                }
            }

           // let us replace all the neighbours
           for (auto neighbour : all_neighbours)
           {
               const auto neigh_side = neighbour.Side();
               auto neigh_el = neighbour.Element();

               /* 
               let us take into account the possibility that
               one triangle might be neighbour of two cylinders
               */

               if (!neigh_el->IsGeoBlendEl())
               {
                   const auto neigh_index = neigh_el->Index();
                   TPZChangeEl::ChangeToGeoBlend(&gmesh, neigh_index);
               }
               else
               {
                   neigh_el->SetNeighbourForBlending(neigh_side);
               }
           }
        }
    }

    for (auto cylinder : found_cylinders)
    {
        if (!cylinder.second)
        {
            PZError<<__PRETTY_FUNCTION__
            <<"\n cylinder " << cylinder.first << " not found in mesh" << std::endl;
        }
    }
}

void TPZMeshOperator::printVTKWJacInfo(std::string filename, TPZGeoMesh* gmesh) {
    TPZVec<REAL> elData(gmesh->NElements(), -100);
    for (int i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl* gel = gmesh->Element(i);
        if(!gel) DebugStop();
        const int geldim = gel->Dimension();
        TPZManVector<REAL,3> qsi(geldim,0.);
        gel->CenterPoint(gel->NSides()-1, qsi);
        REAL detjac = -1000;
        TPZFMatrix<REAL> jac(3,3,0.), axes(3,3,0.), jacinv(3,3,0.);
        gel->Jacobian(qsi, jac, axes, detjac, jacinv);
        elData[i] = detjac;
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, elData);
}
