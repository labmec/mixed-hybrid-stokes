#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fBcNormalVec.clear(); fBcNormalVec.reserve(10);
    fBcTangentialVec.clear(); fBcTangentialVec.reserve(10);
    fBcAxisymmetryVec.clear(); fBcAxisymmetryVec.reserve(10);
    fDomain.clear(); fDomain.reserve(10);
    fAxisymmetryDomain.clear(); fAxisymmetryDomain.reserve(10);
    farcData.clear();
    fcylData.clear();
}

// deconstructor
ProblemData::~ProblemData(){
    
}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file){
    std::ifstream filejson(file);

    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    // checking infos in the json file
    if(input.find("MeshName") == input.end()) DebugStop();
    if(input.find("CreateMsh") == input.end()) DebugStop();
    if(input.find("HdivType")==input.end()) DebugStop();
    if(input.find("isAxisymmetric")==input.end()) DebugStop();
    if(input.find("VelpOrder") == input.end()) DebugStop();
    if(input.find("TracpOrder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Resolution")==input.end()) DebugStop();
    if(input.find("StaticCondensation") == input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
    if(input.find("NormalBoundary") == input.end()) DebugStop();
    if(input.find("TangentialBoundary") == input.end()) DebugStop();
    if(input.find("AxisymmetryDomain") == input.end()) DebugStop();
    if(input.find("AxisymmetryBoundary") == input.end()) DebugStop();
    if(input.find("InterfaceID") == input.end()) DebugStop();
    if(input.find("LambdaID")==input.end()) DebugStop();
    if(input.find("AxiLambdaID") == input.end()) DebugStop();
    if(input.find("AxiInterfaceID")==input.end()) DebugStop();
    if(input.find("FluxInterfaceID") == input.end()) DebugStop();
    if(input.find("HasAnalyticSolution") == input.end()) DebugStop();
        
    // accessing and assigning values
    fMeshName = input["MeshName"];
    
    fMshFile = input["CreateMsh"];
    
    fHdivtype = input["HdivType"];
    
    fAxisymmetric = input["isAxisymmetric"];

    fVelpOrder = input["VelpOrder"];
    
    fTracpOrder = input["TracpOrder"];
    
    fDim = input["Dim"];
    
    fResolution = input["Resolution"];
    
    fCondensedElement = input["StaticCondensation"];
    
    fhasAnalyticSolution = input["HasAnalyticSolution"];
    
    if (input.find("ObstructionID") != input.end())
        fObstructionID = input["ObstructionID"];
    
    if (input.find("CsvFile") != input.end())
        fcsvFile = input["CsvFile"];

    DomainData domaindata;
    for(auto& domainjson : input["Domain"]){
        if(domainjson.find("name") == domainjson.end()) DebugStop();
        if(domainjson.find("matID") == domainjson.end()) DebugStop();
        if(domainjson.find("viscosity") == domainjson.end()) DebugStop();
        
        domaindata.name = domainjson["name"];
        domaindata.matID = domainjson["matID"];
        domaindata.viscosity = domainjson["viscosity"];
        
        fDomain.push_back(domaindata);
    }
    
    BcData bcNormaldata;
    for(auto& bcjson : input["NormalBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcNormaldata.name = bcjson["name"];
        bcNormaldata.type = bcjson["type"];
        bcNormaldata.matID = bcjson["matID"];
        bcNormaldata.value[0] = bcjson["value"];
        
        fBcNormalVec.push_back(bcNormaldata);
    }
    
    BcData bcTangentialdata;
    for(auto& bcjson : input["TangentialBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcTangentialdata.name = bcjson["name"];
        bcTangentialdata.type = bcjson["type"];
        for (int i = 0; i < fDim; i++)
            bcTangentialdata.value[i] = bcjson["value"][i];
        bcTangentialdata.matID = bcjson["matID"];
        
        fBcTangentialVec.push_back(bcTangentialdata);
    }

    if (fAxisymmetric)
    {
        AxisymmetryDomainData axidomaindata;
        for(auto& axidomainjson : input["AxisymmetryDomain"])
        {
            if(axidomainjson.find("name") == axidomainjson.end()) DebugStop();
            if(axidomainjson.find("matID") == axidomainjson.end()) DebugStop();
            if(axidomainjson.find("viscosity") == axidomainjson.end()) DebugStop();
            if(axidomainjson.find("radius") == axidomainjson.end()) DebugStop();
            
            axidomaindata.name = axidomainjson["name"];
            axidomaindata.matID = axidomainjson["matID"];
            axidomaindata.viscosity = axidomainjson["viscosity"];
            axidomaindata.radius = axidomainjson["radius"];
            
            fAxisymmetryDomain.push_back(axidomaindata);
        }

        BcData bcAxisymmetrydata;
        for(auto& bcjson : input["AxisymmetryBoundary"])
        {
            if(bcjson.find("name") == bcjson.end()) DebugStop();
            if(bcjson.find("type") == bcjson.end()) DebugStop();
            if(bcjson.find("value") == bcjson.end()) DebugStop();
            if(bcjson.find("matID") == bcjson.end()) DebugStop();

            bcAxisymmetrydata.name = bcjson["name"];
            bcAxisymmetrydata.type = bcjson["type"];
            bcAxisymmetrydata.value[0] = bcjson["value"];
            bcAxisymmetrydata.matID = bcjson["matID"];
            
            fBcAxisymmetryVec.push_back(bcAxisymmetrydata);
        }
    }
    
    fInterfaceID = input["InterfaceID"];
    fLambdaID = input["LambdaID"];
    fAxiLambdaID = input["AxiLambdaID"];
    fAxiInterfaceID = input["AxiInterfaceID"];
    fFluxInterfaceID = input["FluxInterfaceID"];

    if (fCondensedElement && fHdivtype != EConstant) fMeshVector.resize(4);
    else fMeshVector.resize(2);
}

void ProblemData::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh Name: " << fMeshName << std::endl << std::endl;
    
    out << "Hdiv Type: " << fHdivtype << std::endl << std::endl;

    out << "IsAxisymmetric: " << fAxisymmetric << std::endl << std::endl;
    
    out << "Velocity pOrder: " << fVelpOrder << std::endl << std::endl;
    
    out << "Traction pOrder: " << fTracpOrder << std::endl << std::endl;
    
    out << "Dimension: " << fDim << std::endl << std::endl;
    
    out << "Resolution: " << fResolution << std::endl << std::endl;
    
    out << "Static Condensation: " << fCondensedElement << std::endl << std::endl;
    
    out << "Domain: " << std::endl;
    
    for(const auto& domaindata : fDomain){
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain Viscosity: " << domaindata.viscosity << std::endl << std::endl;
    }
    
    out << "Normal Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fBcNormalVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Tangential Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fBcTangentialVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }

    out << "Axisymmetric Domain at r=0: " << std::endl;
    
    for(const auto& domaindata : fAxisymmetryDomain){
        out << "  Axisymmetric Domain Name: " << domaindata.name << std::endl;
        out << "  Axisymmetric Domain MatID: " << domaindata.matID << std::endl;
        out << "  Axisymmetric Domain Viscosity: " << domaindata.viscosity << std::endl << std::endl;
        out << "  Axisymmetric Domain radius: " << domaindata.radius << std::endl << std::endl;
    }

    out << "Axisymmetric Boundary Conditions at r=0: " << std::endl;
    
    for(const auto& bcdata : fBcAxisymmetryVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }

    out << "Interface Elements ID: " << fInterfaceID << std::endl << std::endl;
    
    out << "Lambda Elements ID: " << fLambdaID << std::endl << std::endl;
}

void ProblemData::ReadCirclesData()
{
    std::ifstream file;
    file.open(fcsvFile);
    
    std::string line;
    
    while (std::getline(file, line))
    {
        ArcData arc;
        std::string tempString;
        
        std::stringstream inputData(line);
        
        std::getline(inputData, tempString, ',');
        arc.radius = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        arc.xc = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        arc.yc = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        arc.zc  = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        arc.matID = atoi(tempString.c_str());
        
        farcData.push_back(arc);
    }
}

void ProblemData::ReadCylindersData()
{
    std::ifstream file;
    file.open(fcsvFile);
    
    std::string line;
    
    while (std::getline(file, line))
    {
        CylinderData cyl;
        std::string tempString;
        
        std::stringstream inputData(line);
        
        std::getline(inputData, tempString, ',');
        cyl.radius = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.xc = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.yc = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.zc = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.xaxis = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.yaxis = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.zaxis = atof(tempString.c_str());
        
        std::getline(inputData, tempString, ',');
        cyl.matID = atoi(tempString.c_str());

        fcylData.push_back(cyl);
    }
}
