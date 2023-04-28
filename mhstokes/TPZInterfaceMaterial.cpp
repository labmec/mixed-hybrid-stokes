#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzlog.h>

#include "TPZInterfaceMaterial.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesInterface");
#endif


TPZInterfaceMaterial::TPZInterfaceMaterial(int matID, int dimension) : TBase(matID), fdimension(dimension) {
    
}

TPZInterfaceMaterial::~TPZInterfaceMaterial(){
    
}

STATE TPZInterfaceMaterial::InnerProductVec(TPZFMatrix<STATE>& S, TPZFMatrix<STATE>& T){
#ifdef DEBUG
    if(S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Cols()) DebugStop();
#endif
    
    STATE val = 0.;
    
    for(int j=0; j<S.Cols(); j++){
        for(int i=0; i<S.Rows(); i++){
            val += S(i, j)*T(i, j);
        }
    }
    
    return val;
}

void TPZInterfaceMaterial::ContributeInterface(const TPZMaterialDataT<STATE>& data, const std::map<int, TPZMaterialDataT<STATE>>& dataleft, const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) {
    
    if(dataleft.find(fVindex) == dataleft.end()) DebugStop();
    if(dataright.find(fPindex) == dataright.end()) DebugStop();
    
    const TPZMaterialDataT<STATE>& vDataLeft = dataleft.find(fVindex)->second;
    const TPZMaterialDataT<STATE>& pDataRight = dataright.find(fPindex)->second;
    
    const TPZFNMatrix<9, REAL>& tan = pDataRight.axes;
    
    int64_t nShapeV = vDataLeft.fVecShapeIndex.NElements();
    int64_t nShapeLambda = pDataRight.phi.Rows();
    
    int nStateVariablesL = fdimension - 1;
    
    int64_t VRows = vDataLeft.fDeformedDirections.Rows();
    int64_t VCols = vDataLeft.fDeformedDirections.Cols();
    
    TPZFNMatrix<3, REAL> phiV(VRows, VCols, 0.);
    
    if(vDataLeft.fNeedsDeformedDirectionsFad){
        for(int row=0; row<VRows; row++){
            for(int col=0; col<VCols; col++){
                phiV(row, col) = vDataLeft.fDeformedDirectionsFad(row, col).val();
            }
        }
    } else {
        DebugStop();
    }
    
    for(int vFunction_i=0; vFunction_i<nShapeV; vFunction_i++){
        STATE LambdaDotPhiV = 0.;
        
        TPZFNMatrix<3, STATE> phiVi(3,1,0.);
        for(int row=0; row<3; row++){
            phiVi(row, 0) = phiV(row, vFunction_i);
        }
        
        for(int LambdaFunction_j=0; LambdaFunction_j<nShapeLambda; LambdaFunction_j++){
            TPZFNMatrix<3, STATE> lambda_j(3,1,0.);
            TPZFNMatrix<3, STATE> phiLambda = pDataRight.phi;
            
            for(int f=0; f<nStateVariablesL; f++){
                lambda_j.Zero();
                for(int row=0; row<3; row++){
                    lambda_j(row, 0) += phiLambda(LambdaFunction_j,0)*tan(f, row);
                }
                
                STATE fact = fMultiplier*InnerProductVec(phiVi, lambda_j)*weight;
                ek(vFunction_i, LambdaFunction_j*nStateVariablesL+f+nShapeV) += fact;
                ek(LambdaFunction_j*nStateVariablesL+f+nShapeV, vFunction_i) += fact;
            }
        }
    }


#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        std::stringstream sout;
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZInterfaceMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    
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
