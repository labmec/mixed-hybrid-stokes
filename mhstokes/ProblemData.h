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
    
    bool fMshFile = false;
    
    std::string fHdivtype = "none";
    
    int fVelpOrder = -1; // polynomial approximation order for velocity
    
    int fTracpOrder = -1; // polynomial approximation order for traction
    
    int fDim = -1;
    
    int fresolution = -1;
    
    bool fCondensedElement = false;
    
    // struct responsible to summarize all the data from every domain
    struct fDomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        double viscosity = -1.; // domain viscosity
    };
    
    std::vector<fDomainData> fdomain; // vector containing every domain created
    
    // struct responsible to summarize all the data from every velocity boundary condition
    struct fBcNormalData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<double,1>  value = {0}; // bc value
        int matID = 0; // bc material ID
    };

    // struct responsible to summarize all the data from every traction boundary condition
    struct fBcTangentialData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<double,1> value = {0}; // bc value
        int matID = 0; // bc material ID
    };
    
    std::vector<fBcNormalData> fbcNormalvec; // vector containg all the velocity bcs info
    
    std::vector<fBcTangentialData> fbcTangentialvec; // vector containg all the traction bcs info
    
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
    
    const bool& CreateMshFile() const {return fMshFile;}
    bool& CreateMshFile() {return fMshFile;}
    
    const std::string& HdivType() const {return fHdivtype;}
    std::string& HdivType() {return fHdivtype;}
    
    const int& VelpOrder() const {return fVelpOrder;}
    int& VelpOrder(){return fVelpOrder;}
    
    const int& TracpOrder() const {return fTracpOrder;}
    int& tracpOrder(){return fTracpOrder;}
    
    const int& Dim() const {return fDim;}
    int& Dim(){return fDim;}
    
    const int& Resolution() const {return fresolution;}
    int& Resolution() {return fresolution;}
    
    const bool& CondensedElements() const {return fCondensedElement;}
    bool& CondensedElements() {return fCondensedElement;}

    const std::vector<fDomainData>& DomainVec() const {return fdomain;}
    std::vector<fDomainData>& DomainVec(){return fdomain;}
    
    const std::vector<fBcNormalData>& NormalBCs() const {return fbcNormalvec;}
    std::vector<fBcNormalData>& NormalBCs(){return fbcNormalvec;}
    
    const std::vector<fBcTangentialData>& TangentialBCs() const {return fbcTangentialvec;}
    std::vector<fBcTangentialData>& TangentialBCs(){return fbcTangentialvec;}
    
    const int& InterfaceID() const{return finterfaceID;}
    int& InterfaceID(){return finterfaceID;}
    
    const int& LambdaID() const{return flambdaID;}
    int& LambdaID(){return flambdaID;}
    
    const TPZVec<TPZCompMesh*>& MeshVector() const {return fMeshVectorr;}
    TPZVec<TPZCompMesh*>& MeshVector() {return fMeshVectorr;}
};

#endif
