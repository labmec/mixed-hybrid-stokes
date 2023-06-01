#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fbcNormalvec.clear();
    fbcTangentialvec.clear();
    fdomain.clear();
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
    if(input.find("VelpOrder") == input.end()) DebugStop();
    if(input.find("TracpOrder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Resolution")==input.end()) DebugStop();
    if(input.find("StaticCondensation") == input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
    if(input.find("NormalBoundary") == input.end()) DebugStop();
    if(input.find("TangentialBoundary") == input.end()) DebugStop();
    if(input.find("InterfaceID") == input.end()) DebugStop();
    if(input.find("LambdaID")==input.end()) DebugStop();
    
        
    // accessing and assigning values
    fMeshName = input["MeshName"];
    
    fMshFile = input["CreateMsh"];
    
    fHdivtype = input["HdivType"];
    
    fVelpOrder = input["VelpOrder"];
    
    fTracpOrder = input["TracpOrder"];
    
    fDim = input["Dim"];
    
    fresolution = input["Resolution"];
    
    fCondensedElement = input["StaticCondensation"];
    
    fDomainData domaindata;
    for(auto& domainjson : input["Domain"]){
        if(domainjson.find("name") == domainjson.end()) DebugStop();
        if(domainjson.find("matID") == domainjson.end()) DebugStop();
        if(domainjson.find("viscosity") == domainjson.end()) DebugStop();
        
        domaindata.name = domainjson["name"];
        domaindata.matID = domainjson["matID"];
        domaindata.viscosity = domainjson["viscosity"];
        
        fdomain.push_back(domaindata);
    }
    
    fBcNormalData bcNormaldata;
    for(auto& bcjson : input["NormalBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcNormaldata.name = bcjson["name"];
        bcNormaldata.type = bcjson["type"];
        bcNormaldata.matID = bcjson["matID"];
        bcNormaldata.value[0] = bcjson["value"];
        
        fbcNormalvec.push_back(bcNormaldata);
    }
    
    fBcTangentialData bcTangentialdata;
    for(auto& bcjson : input["TangentialBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcTangentialdata.name = bcjson["name"];
        bcTangentialdata.type = bcjson["type"];
        bcTangentialdata.value[0] = bcjson["value"];
        bcTangentialdata.matID = bcjson["matID"];
        
        fbcTangentialvec.push_back(bcTangentialdata);
    }
    
    finterfaceID = input["InterfaceID"];
    flambdaID = input["LambdaID"];
    
}

void ProblemData::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh Name: " << fMeshName << std::endl << std::endl;
    
    out << "Hdiv Type: " << fHdivtype << std::endl << std::endl;
    
    out << "Velocity pOrder: " << fVelpOrder << std::endl << std::endl;
    
    out << "Traction pOrder: " << fTracpOrder << std::endl << std::endl;
    
    out << "Dimension: " << fDim << std::endl << std::endl;
    
    out << "Resolution: " << fresolution << std::endl << std::endl;
    
    out << "Static Condensation: " << fCondensedElement << std::endl << std::endl;
    
    out << "Domain: " << std::endl;
    
    for(const auto& domaindata : fdomain){
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain Viscosity: " << domaindata.viscosity << std::endl << std::endl;
    }
    
    out << "Normal Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fbcNormalvec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Tangential Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fbcTangentialvec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Interface Elements ID: " << finterfaceID << std::endl << std::endl;
    
    out << "Lambda Elements ID: " << flambdaID << std::endl << std::endl;
}

