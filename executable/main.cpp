#include <iostream>
#include <fstream>
#include <string>
#include <json.hpp>

#include"ProblemData.h"

int main(){
    std::string filenamejson =  "StokesData.json";
    
    ProblemData simData;

    simData.ReadJson(filenamejson);

    simData.PrintData();
    
    // ------------------------ Getting number of domains and fractures ------------------------
//    if(input.find("Domains") == input.end()) DebugStop();
//    const int ndom = input["Domains"].size();
//    if(input.find("Fractures") == input.end());
//    const int nfrac = input["Fractures"].size();
//    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(ndom+nfrac+1);
//    int countPhi = 0;
    
	return 0;
}
