#include <iostream>
#include <string>
#include<fstream>
#include <json.hpp>
#include <pzerror.h>

#include "ProblemData.h"

using json = nlohmann::json;
using namespace std;

ProblemData::ProblemData(){
    fbcvec.clear();
}

void ProblemData::ReadJson(std::string file){
//        std::ifstream filejson(file);
//        json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
//
//        if(input.find("TolDist") == input.end()) std::cout << "Ferrou!" << std::endl;
//        double toldist = input["TolDist"];
//        std::cout << "toldist = " << toldist << std::endl;
//
//        json domains = input["Domains"];
//        std::cout << "domains = " << domains << std::endl;
//
    std::ifstream filejson(file);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
        
    for(auto& bcjson : input["Boundary"]){
        const string name = bcjson["name"];
        const int type = bcjson["type"];
        const double value = bcjson["value"];
        const int matid = bcjson["matID"];
        
        fbcvec.push_back(std::make_tuple(name,type,value,matid));
    }
    
    for(auto& mytup : fbcvec){
        cout << std::get<0>(mytup) << endl;
        cout << std::get<1>(mytup) << endl;
        cout << std::get<2>(mytup) << endl;
        cout << std::get<3>(mytup) << endl;
    }
    DebugStop();
    

//    Mesh = Input["Mesh"];
//    nOrder = Input["nOrder"];
//    Domain = Input["Domain"];
//    Boundary = Input["Boundary"];
}

void ProblemData::PrintData(){
    std::cout << "Simulation Mesh " << fMeshName << std::endl;
    std::cout << "Simulation nOrder " << fpOrder << std::endl;

}
