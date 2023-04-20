#include <pzfmatrix.h>
#include <TPZBndCondT.h>

#include "TPZInterfaceMaterial.h"

TPZInterfaceMaterial::TPZInterfaceMaterial(int matID, int dimension) : TBase(matID), fdimension(dimension) {
    
}

TPZInterfaceMaterial::~TPZInterfaceMaterial(){
    
}

void TPZInterfaceMaterial::ContributeInterface(const TPZMaterialDataT<STATE>& data, const std::map<int, TPZMaterialDataT<STATE>>& dataleft, const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) {
    
    DebugStop();
        
}

void TPZInterfaceMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &dataleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    
    DebugStop();
}

void TPZInterfaceMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TPZInterfaceMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  REAL weight, TPZFMatrix<STATE> &ek,
                  TPZFMatrix<STATE> &ef,
                  TPZBndCondT<STATE> &bc){
    DebugStop();
}

void TPZInterfaceMaterial::FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, std::map<int, TPZMaterialDataT<STATE>> &datavec_right) {
    
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsSol = true;
    datavec_right[1].fNeedsNormal = true;
    datavec_right[1].fNeedsSol = true;
    datavec_left[0].fNeedsDeformedDirectionsFad = true;
}

void TPZInterfaceMaterial::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                       const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                       const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                       int var, TPZVec<STATE> &Solout) {
    DebugStop();
}


void TPZInterfaceMaterial::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                       const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                       const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                       int var, TPZVec<STATE> &Solout,
                       TPZCompEl *left,TPZCompEl *right) {
    DebugStop();
}

void TPZInterfaceMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol) {
    
    DebugStop();
}
