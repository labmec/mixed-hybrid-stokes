#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fbcvec.clear();
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
    if(input.find("Boundary") == input.end()) DebugStop();
    if(input.find("InterfaceID") == input.end()) DebugStop();
    
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
    
    fBcData bcdata;
    for(auto& bcjson : input["Boundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcdata.name = bcjson["name"];
        bcdata.type = bcjson["type"];
        bcdata.value = bcjson["value"];
        bcdata.matID = bcjson["matID"];
        
        fbcvec.push_back(bcdata);
    }
    
    finterfaceID = input["InterfaceID"];
    
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
    
    out << "Boundary COnditions: " << std::endl;
    
    for(const auto& bcdata : fbcvec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Interface Elements ID: " << finterfaceID << std::endl << std::endl;
}

