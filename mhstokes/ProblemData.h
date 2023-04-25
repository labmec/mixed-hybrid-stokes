#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include <json.hpp>
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData{
private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";
    
    int fVelpOrder = -1; // polynomial approximation order for velocity
    
    int fTracpOrder = -1; // polynomial approximation order for traction
    
    int fDim = -1;
    
    // struct responsible to summarize all the data from every domain
    struct fDomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        double viscosity = -1.; // domain viscosity
    };
    
    std::vector<fDomainData> fdomain; // vector containing every domain created
    
    // struct responsible to summarize all the data from every velocity boundary condition
    struct fBcVelData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<double,3>  value = {0,0,0}; // bc value
        int matID = 0; // bc material ID
    };

    // struct responsible to summarize all the data from every traction boundary condition
    struct fBcTracData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<double,3> value = {0,0,0}; // bc value
        int matID = 0; // bc material ID
    };
    
    std::vector<fBcVelData> fbcvelvec; // vector containg all the velocity bcs info
    
    std::vector<fBcTracData> fbctracvec; // vector containg all the traction bcs info
    
    int finterfaceID = -1;
    
    int flambdaID = -1;
    
    TPZVec<TPZCompMesh*> fMeshVectorr;
    
public:
    ProblemData();
    
    ~ProblemData();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout);
    // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;} //setter using reference variable;
    std::string& MeshName(){return fMeshName;}  //getter using reference variable
    
    const int& VelpOrder() const {return fVelpOrder;}
    int& VelpOrder(){return fVelpOrder;}
    
    const int& TracpOrder() const {return fTracpOrder;}
    int& tracpOrder(){return fTracpOrder;}
    
    const int& Dim() const {return fDim;}
    int& Dim(){return fDim;}

    const std::vector<fDomainData>& DomainVec() const {return fdomain;}
    std::vector<fDomainData>& DomainVec(){return fdomain;}
    
    const std::vector<fBcVelData>& VelBCs() const {return fbcvelvec;}
    std::vector<fBcVelData>& VelBCs(){return fbcvelvec;}
    
    const std::vector<fBcTracData>& TracBCs() const {return fbctracvec;}
    std::vector<fBcTracData>& TracBCs(){return fbctracvec;}
    
    const int& InterfaceID() const{return finterfaceID;}
    int& InterfaceID(){return finterfaceID;}
    
    const int& LambdaID() const{return flambdaID;}
    int& LambdaID(){return flambdaID;}
    
    const TPZVec<TPZCompMesh*>& MeshVector() const {return fMeshVectorr;}
    TPZVec<TPZCompMesh*>& MeshVector() {return fMeshVectorr;}
};

#endif
