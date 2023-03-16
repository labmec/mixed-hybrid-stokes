#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
SimulationData::SimulationData(){
    fbcvec.clear();
    fdomain.clear();
}

// deconstructor
SimulationData::~SimulationData(){
    
}

// readjson function. takes a json function as parameter and completes the required simulation data
void SimulationData::ReadJson(std::string file){
    std::ifstream filejson(file);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    // checking infos in the json file
    if(input.find("MeshName") == input.end()) DebugStop();
    if(input.find("pOrder") == input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
    if(input.find("Boundary") == input.end()) DebugStop();
    
    // accessing and assigning values
    fMeshName = input["MeshName"];
    
    fpOrder = input["pOrder"];
    
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
    
}

void SimulationData::Print(std::ostream& out){
    out << "A new simulation has been started: \n\n";
    out << "Simulation Mesh Name: " << fMeshName << std::endl << std::endl;
    
    out << "Simulation pOrder: " << fpOrder << std::endl << std::endl;
    
    out << "Simulation Domain: " << std::endl;
    
    for(auto domaindata : fdomain){
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain Viscosity: " << domaindata.viscosity << std::endl << std::endl;
    }
    
    out << "Simulation Boundary COnditions: " << std::endl;
    
    for(auto bcdata : fbcvec){
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
}
