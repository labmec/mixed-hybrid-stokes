#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzlog.h>

#include "TPZInterfaceAxisymStokesMaterial.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesInterface");
#endif


TPZInterfaceAxisymStokesMaterial::TPZInterfaceAxisymStokesMaterial(int matID, int dimension, bool isAxisymmetric) : TBase(matID), fdimension(dimension), faxisymmetry(isAxisymmetric) {
    
}

TPZInterfaceAxisymStokesMaterial::~TPZInterfaceAxisymStokesMaterial(){
    
}

STATE TPZInterfaceAxisymStokesMaterial::InnerProductVec(TPZFMatrix<STATE>& S, TPZFMatrix<STATE>& T){
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

void TPZInterfaceAxisymStokesMaterial::ContributeInterface(const TPZMaterialDataT<STATE>& data, const std::map<int, TPZMaterialDataT<STATE>>& dataleft, const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) {
    
    if(dataleft.find(fVindex) == dataleft.end()) DebugStop();
    if(dataright.find(fPindex) == dataright.end()) DebugStop();
    
    const TPZMaterialDataT<STATE>& vDataLeft = dataleft.find(fVindex)->second;
    const TPZMaterialDataT<STATE>& pDataRight = dataright.find(fPindex)->second;
    
    const TPZFNMatrix<9, REAL>& axes = pDataRight.axes;
    
    int64_t nShapeV = vDataLeft.fVecShapeIndex.NElements(); // number of Hdiv velocity/displacement shape functions
    int64_t nShapeLambda = pDataRight.phi.Rows(); // number of lambda shape functions
    
    TPZFNMatrix<100, REAL> PhiV(3, nShapeV, 0.0);
    TPZFNMatrix<100, REAL> PhiLambdaT(3, fdimension * nShapeLambda, 0.0); //lambda in the tangential directions
    
    if (vDataLeft.fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {
            for (int64_t k = 0; k < 3; k++)
            {
                PhiV(k, j) = vDataLeft.fDeformedDirectionsFad(k, j).val();
            }
        }
    }

    for (int64_t j = 0; j < nShapeLambda; j++)
    {
        for (int k = 0; k < fdimension; k++)
        {
            TPZManVector<REAL,3> tangent(3,0.0);
            for (int i = 0; i < 3; i++) tangent[i] = axes(k,i);
            for (int64_t i = 0; i < 3; i++)
            {
                
                PhiLambdaT(i, fdimension*j+k) = pDataRight.phi(j, 0) * axes(k, i); // lambda shape function at tangential direction
            }
        }
    }

    REAL factor = fMultiplier * weight;
    ek.AddContribution(0, nShapeV, PhiV, true, PhiLambdaT, false, factor);
    ek.AddContribution(nShapeV, 0, PhiLambdaT, true, PhiV, false, factor);

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

void TPZInterfaceAxisymStokesMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    
    DebugStop();
    
}

void TPZInterfaceAxisymStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    DebugStop();
}

void TPZInterfaceAxisymStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  REAL weight, TPZFMatrix<STATE> &ek,
                  TPZFMatrix<STATE> &ef,
                  TPZBndCondT<STATE> &bc){
    DebugStop();
}

void TPZInterfaceAxisymStokesMaterial::FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, std::map<int, TPZMaterialDataT<STATE>> &datavec_right) {
    
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsSol = true;
    datavec_right[1].fNeedsNormal = true;
    datavec_right[1].fNeedsSol = true;
    datavec_left[0].fNeedsDeformedDirectionsFad = true;
}

void TPZInterfaceAxisymStokesMaterial::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                       const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                       const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                       int var, TPZVec<STATE> &Solout) {
    DebugStop();
}


void TPZInterfaceAxisymStokesMaterial::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                       const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                       const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                       int var, TPZVec<STATE> &Solout,
                       TPZCompEl *left,TPZCompEl *right) {
    DebugStop();
}

void TPZInterfaceAxisymStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol) {
    
    DebugStop();
}

int TPZInterfaceAxisymStokesMaterial::GetIntegrationOrder(const TPZVec<int> &porder_left, const TPZVec<int> &porder_right) const
{
    int maxl = 0, maxr = 0;
    for (auto porder: porder_left) {
        maxl = std::max(maxl,porder);
    }
    for (auto porder: porder_right) {
        maxr = std::max(maxr,porder);
    }
    return (maxl+maxr);
}
