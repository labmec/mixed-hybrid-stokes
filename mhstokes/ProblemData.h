#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class SimulationData{
private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none"; // name of the finite element geometric mesh
    int fpOrder = -1; // polynomial approximation order
    
    // struct responsible to summarize all the data from every domain
    struct fDomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        double viscosity = -1.; // domain viscosity
    };
    
    std::vector<fDomainData> fdomain; // vector containing every domain created
    
    // struct responsible to summarize all the data from every single boundary condition
    struct fBcData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        double value = 0.; // bc value
        int matID = 0; // bc material ID
    };
    
    std::vector<fBcData> fbcvec; // vector containg all the bcs info
    
public:
    SimulationData();
    
    ~SimulationData();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout); // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;} // setter using reference variable
    std::string& MeshName(){return fMeshName;} // getter using reference variable
    
    const int& POrder() const {return fpOrder;}
    int& POrder(){return fpOrder;}

    const std::vector<fDomainData>& DomainVec() const {return fdomain;}
    std::vector<fDomainData>& DomainVec(){return fdomain;}
    
    const std::vector<fBcData>& BCs() const {return fbcvec;}
    std::vector<fBcData>& BCs(){return fbcvec;}
    
};

#endif
