#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include <json.hpp>
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
    // struct responsible to summarize all the data from every domain
    struct DomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        double viscosity = -1.; // domain viscosity
    };
    
    // struct responsible to summarize all the data from axisymmetric infinitesimal tube at r=0
    struct AxisymmetryDomainData {
        std::string name = "none"; // domains name
        int matID = 0; // bc material ID
        double viscosity = -1.; // domain viscosity
        REAL radius = 0; // domain radius
    };

    // struct responsible to store boundary condition data
    struct BcData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<double,1>  value = {0}; // bc value
        int matID = 0; // bc material ID
    };

private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";
    
    bool fMshFile = false;
    
    int fHdivtype = 0; // 0 -> Standard & 1 -> Constant
    
    int fVelpOrder = -1; // polynomial approximation order for velocity
    
    int fTracpOrder = -1; // polynomial approximation order for traction
    
    int fDim = -1;
    
    int fResolution = -1;

    int fAxisymmetric = 0; // whether it is an axisymmetric or cartesian simulation. 0 - cartesian, 1 - axisymmetric
    
    bool fCondensedElement = false;
    
    std::vector<DomainData> fDomain; // vector containing every domain created

    std::vector<AxisymmetryDomainData> fAxisymmetryDomain; // vector containing every axisymmetric domain created at r=0
    
    std::vector<BcData> fBcNormalVec; // vector containg all the velocity bcs info
    
    std::vector<BcData> fBcTangentialVec; // vector containg all the traction bcs info

    std::vector<BcData> fBcAxisymmetryVec; //vector containing all axisymmetry boundary info at r=0
    
    int fInterfaceID = -1;
    
    int fLambdaID = -1;

    int fAxiLambdaID = -1;

    int fAxiInterfaceID = -1;

    int fFluxInterfaceID = -1;
    
    TPZVec<TPZCompMesh*> fMeshVector;
    
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
    
    const int& HdivType() const {return fHdivtype;}
    int& HdivType() {return fHdivtype;}

    const int& Axisymmetric() const {return fAxisymmetric;}
    int& Axisymmetric() {return fAxisymmetric;}
    
    const int& VelpOrder() const {return fVelpOrder;}
    int& VelpOrder(){return fVelpOrder;}
    
    const int& TracpOrder() const {return fTracpOrder;}
    int& tracpOrder(){return fTracpOrder;}
    
    const int& Dim() const {return fDim;}
    int& Dim(){return fDim;}
    
    const int& Resolution() const {return fResolution;}
    int& Resolution() {return fResolution;}
    
    const bool& CondensedElements() const {return fCondensedElement;}
    bool& CondensedElements() {return fCondensedElement;}

    const std::vector<DomainData>& DomainVec() const {return fDomain;}
    std::vector<DomainData>& DomainVec(){return fDomain;}

    const std::vector<AxisymmetryDomainData>& AxisymmetryDomainVec() const {return fAxisymmetryDomain;}
    std::vector<AxisymmetryDomainData>& AxisymmetryDomainVec(){return fAxisymmetryDomain;}
    
    const std::vector<BcData>& NormalBCs() const {return fBcNormalVec;}
    std::vector<BcData>& NormalBCs(){return fBcNormalVec;}
    
    const std::vector<BcData>& TangentialBCs() const {return fBcTangentialVec;}
    std::vector<BcData>& TangentialBCs(){return fBcTangentialVec;}

    const std::vector<BcData>& AxisymmetryBCs() const {return fBcAxisymmetryVec;}
    std::vector<BcData>& AxisymmetryBCs(){return fBcAxisymmetryVec;}
    
    const int& InterfaceID() const{return fInterfaceID;}
    int& InterfaceID(){return fInterfaceID;}
    
    const int& LambdaID() const{return fLambdaID;}
    int& LambdaID(){return fLambdaID;}

    const int& AxiLambdaID() const{return fAxiLambdaID;}
    int& AxiLambdaID(){return fAxiLambdaID;}
    
    const int& AxiInterfaceID() const{return fAxiInterfaceID;}
    int& AxiInterfaceID(){return fAxiInterfaceID;}

    const int& FluxInterfaceID() const{return fFluxInterfaceID;}
    int& FluxInterfaceID(){return fFluxInterfaceID;}
    
    const TPZVec<TPZCompMesh*>& MeshVector() const {return fMeshVector;}
    TPZVec<TPZCompMesh*>& MeshVector() {return fMeshVector;}
};

#endif
