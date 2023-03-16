#include <iostream>
#include <fstream>
#include <string>
#include <json.hpp>

#include"ProblemData.h"

int main(){
    std::string filenamejson =  "StokesData.json";
    std::ofstream simulationfile("SimulationData.txt");
    
    SimulationData simData;

    simData.ReadJson(filenamejson);
    
    simData.Print(simulationfile);
    
	return 0;
}
