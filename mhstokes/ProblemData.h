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
public: 
    enum HDivType {EStandard, EConstant};
    
private:
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
        TPZManVector<double, 3>  value = {0.0, 0.0, 0.0}; // bc value
        int matID = 0; // bc material ID
    };

private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";
    
    bool fMshFile = false;
    
    HDivType fHdivtype = EStandard; // 0 -> Standard & 1 -> Constant
    
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
    
    int fTanVelID = 15;

    int fAxiLambdaID = -1;

    int fAxiInterfaceID = -1;
    
    int fAxiTanVelID = 25;

    int fFluxInterfaceID = -1;
    
    bool fhasAnalyticSolution = false;
    
    TPZVec<TPZCompMesh*> fMeshVector;
    
    int fObstructionID = -1;
    
public:
    ProblemData();
    
    ~ProblemData();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout);
    // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;} //setter using reference variable;
    void SetMeshName(const std::string& meshname) {fMeshName = meshname;}  //getter using reference variable
    
    const bool& CreateMshFile() const {return fMshFile;}
    void SetCreateMshFile(bool create) {fMshFile = create;}
    
    const HDivType& HdivType() const {return fHdivtype;}
    void SetHdivType(HDivType hdivtype) {fHdivtype = hdivtype;}

    const int& Axisymmetric() const {return fAxisymmetric;}
    void SetAxisymmetric(int axisymmetric) {fAxisymmetric = axisymmetric;}
    
    const int& VelpOrder() const {return fVelpOrder;}
    void SetVelpOrder(int velp) {fVelpOrder = velp;}
    
    const int& TracpOrder() const {return fTracpOrder;}
    void SetTracpOrder(int tracp) {fTracpOrder = tracp;}
    
    const int& Dim() const {return fDim;}
    void SetDim(int dim) {fDim = dim;}
    
    const int& Resolution() const {return fResolution;}
    void SetResolution(int res) {fResolution = res;}
    
    const bool& CondensedElements() const {return fCondensedElement;}
    void SetCondensedElements(bool condense) {fCondensedElement = condense;}

    const std::vector<DomainData>& DomainVec() const {return fDomain;}
    void SetDomainVec(const std::vector<DomainData>& vec) {fDomain = vec;}

    const std::vector<AxisymmetryDomainData>& AxisymmetryDomainVec() const {return fAxisymmetryDomain;}
    void SetAxisymmetryDomainVec(const std::vector<AxisymmetryDomainData>& vec) {fAxisymmetryDomain = vec;}
    
    const std::vector<BcData>& NormalBCs() const {return fBcNormalVec;}
    void SetNormalBCs(const std::vector<BcData>& bcs) {fBcNormalVec = bcs;}
    
    const std::vector<BcData>& TangentialBCs() const {return fBcTangentialVec;}
    void SetTangentialBCs(const std::vector<BcData>& bcs) {fBcTangentialVec = bcs;}

    const std::vector<BcData>& AxisymmetryBCs() const {return fBcAxisymmetryVec;}
    void SetAxisymmetryBCs(const std::vector<BcData>& bcs) {fBcAxisymmetryVec = bcs;}
    
    const int& InterfaceID() const{return fInterfaceID;}
    void SetInterfaceID(int id) {fInterfaceID = id;}
    
    const int& LambdaID() const{return fLambdaID;}
    void SetLambdaID(int id ){fLambdaID = id;}
    
    const int& TanVelID() const{return fTanVelID;}
    void SetTanVelID(int id){fTanVelID = id;}

    const int& AxiLambdaID() const{return fAxiLambdaID;}
    void SetAxiLambdaID(int id) {fAxiLambdaID = id;}
    
    const int& AxiInterfaceID() const{return fAxiInterfaceID;}
    void SetAxiInterfaceID(int id) {fAxiInterfaceID = id;}

    const int& AxiTanVelID() const{return fAxiTanVelID;}
    void SetAxiTanVelID(int id){fAxiTanVelID = id;}

    const int& FluxInterfaceID() const{return fFluxInterfaceID;}
    void SetFluxInterfaceID(int id) {fFluxInterfaceID = id;}
    
    const TPZVec<TPZCompMesh*>& MeshVector() const {return fMeshVector;}
    void SetMeshVector(const TPZVec<TPZCompMesh*>& vec) {fMeshVector = vec;}
    
    const bool& HasAnalyticSolution() const{return fhasAnalyticSolution;}
    void SetHasAnalyticSolution(bool analyticalSol){fhasAnalyticSolution = analyticalSol;}
    
    const int& ObstructionID() const{return fObstructionID;}
    void SetObstruction(int obID){fObstructionID = obID;}
};

#endif
