#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fbcvelvec.clear();
    fbctracvec.clear();
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
    if(input.find("VelpOrder") == input.end()) DebugStop();
    if(input.find("TracpOrder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
    if(input.find("BoundaryVelocity") == input.end()) DebugStop();
    if(input.find("BoundaryTraction") == input.end()) DebugStop();
    if(input.find("InterfaceID") == input.end()) DebugStop();
    if(input.find("LambdaID")==input.end()) DebugStop();
    
    // accessing and assigning values
    fMeshName = input["MeshName"];
    
    fVelpOrder = input["VelpOrder"];
    
    fTracpOrder = input["TracpOrder"];
    
    fDim = input["Dim"];
    
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
    
    fBcVelData bcveldata;
    for(auto& bcjson : input["BoundaryVelocity"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcveldata.name = bcjson["name"];
        bcveldata.type = bcjson["type"];
        bcveldata.value = bcjson["value"];
        bcveldata.matID = bcjson["matID"];
        
        fbcvelvec.push_back(bcveldata);
    }
    
    fBcTracData bctracdata;
    for(auto& bcjson : input["BoundaryTraction"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bctracdata.name = bcjson["name"];
        bctracdata.type = bcjson["type"];
        bctracdata.value = bcjson["value"];
        bctracdata.matID = bcjson["matID"];
        
        fbctracvec.push_back(bctracdata);
    }
    
    finterfaceID = input["InterfaceID"];
    flambdaID = input["LambdaID"];
    
}

void ProblemData::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh address: " << fMeshName << std::endl << std::endl;
    
    out << "Dimension: " << fDim << std::endl << std::endl;
    
    out << "Velocity pOrder: " << fVelpOrder << std::endl << std::endl;
    
    out << "Traction pOrder: " << fTracpOrder << std::endl << std::endl;
    
    out << "Domain: " << std::endl;
    
    for(const auto& domaindata : fdomain){
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain Viscosity: " << domaindata.viscosity << std::endl << std::endl;
    }
    
    out << "Velocity Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fbcvelvec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Traction Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fbctracvec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Interface Elements ID: " << finterfaceID << std::endl << std::endl;
    
    out << "Lambda Elements ID: " << flambdaID << std::endl << std::endl;
}

