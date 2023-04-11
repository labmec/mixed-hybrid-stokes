#include <pzfmatrix.h>
#include <TPZBndCondT.h>

#include "TPZStokesMaterial.h"

TPZStokesMaterial::TPZStokesMaterial(): TBase() {
    
}

TPZStokesMaterial::TPZStokesMaterial(int matID, int dimension, double viscosity) : TBase(matID), fdimension(dimension), fviscosity(viscosity){
    
}

TPZStokesMaterial::~TPZStokesMaterial(){
    
}

template <typename TVar>
TVar TPZStokesMaterial::TensorInnerProduct(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    if(S.Rows() != S.Cols() || T.Rows() != T.Cols() || S.Rows() != T.Rows()) DebugStop();
    
    TVar Val = 0;
    
    for(int i=0; i<S.Cols(); i++){
        for(int j=0; j<S.Rows(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
}

void TPZStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef){
    
    int64_t normVecRows = datavec[fVindex].fDeformedDirections.Rows();
    int64_t normVecCols = datavec[fVindex].fDeformedDirections.Cols();
    
    TPZFNMatrix<3, REAL> NormalVec(normVecRows, normVecCols, 0.);
    
    // check if fvecshapeindex exists
    
    NormalVec = datavec[fVindex].fDeformedDirections;
    
    // Setting the phis (is this the correct variable? or should it be fDeformedDirections?)
    // for V
    TPZFMatrix<REAL>& phiV = datavec[fVindex].phi;
    
    // for P
    TPZFMatrix<REAL>& phiP = datavec[fPindex].phi;
    
    int64_t nShapeV = datavec[fVindex].fVecShapeIndex.NElements();
    int64_t nShapeP = phiP.Rows();
    
    // creating a vector of gradients of the H(div) vectors (why is it that way?)
    TPZManVector<TPZFNMatrix<9, REAL>, 18> GradNormVec(normVecCols);
    for(int i=0; i<normVecRows; i++){
        GradNormVec[i].Redim(3, 3);
    }
    
    // what is this??
    if(datavec[fVindex].fNeedsDeformedDirectionsFad){
        for(int e=0; e<normVecRows; e++){
            for(int s=0; s<normVecCols; s++){
                NormalVec(e,s) = datavec[fVindex].fDeformedDirectionsFad(e,s).val();
            }
        }
        
        TPZFNMatrix<9, REAL> Grad0(3,3,0.);
        for(int s=0; s<normVecCols; s++){
            for(int i=0; i < this->Dimension(); i++){
                for(int j=0; j < this->Dimension(); i++){
                    Grad0(i, j) = datavec[fVindex].fDeformedDirectionsFad(i,s).fastAccessDx(j);
                }
            }
            GradNormVec[s] = Grad0;
        }
    }
    
    // declaring some useful variable. What are they for?
    TPZVec<STATE> Force(3,0.);
    
    TPZFMatrix<STATE> phiVi(3,1,0.);
    TPZFMatrix<STATE> phiVj(3,1,0.);
    
    TPZFNMatrix<100, STATE> divphi;
    TPZFNMatrix<10, STATE> dSolVec = datavec[fVindex].dsol[0];
    
    TPZManVector<STATE, 3> un = datavec[fVindex].sol[0];
    STATE pn = datavec[fPindex].sol[0][0];
    STATE divsol = datavec[fVindex].divsol[0][0];
    
    
    // GetTimeStep ?? Is it necessary? How to do it?
    
    TPZFNMatrix<10, STATE> gradUn = dSolVec;
    
    if(this->HasForcingFunction()){
        TPZFMatrix<STATE> gradU;
        TPZVec<STATE> x(3,0.);
        
        x = datavec[fVindex].x;
        this->ForcingFunction()(x, Force);
    }
    
    // Velocity Depending term
    // phiVi: value of the test function (maybe it isn't anymore)
    // GradVi: gradient of the test function
    for(int i=0; i<nShapeV; i++){
        TPZFNMatrix<9, STATE> GradVi(3,3,0.);
        TPZFNMatrix<9, STATE> DUi(3,3,0.);
        
        for(int row=0; row<3; row++){
            phiVi(row,0) = NormalVec(row,i);
            
            for(int col=0; col<3; col++){
                GradVi(row, col) = GradNormVec[i](row, col);
                GradVi(col, row) = GradVi(row, col);
                
                // symetric gradient of the test function
                DUi(row, col) = 0.5*(GradVi(row, col) + GradVi(col, row));
            }
        }
        
        STATE divUi = datavec[fVindex].divphi(i, 0);
        
        // f source term
        STATE phiDotF = 0.;
        
        for(int row=0; row<3; row++){
            phiDotF += phiVi(row)*Force[row];
        }
        
        ef(i) += weight*phiDotF;
        
        // Stiffness Matrix A (depending term of V)
        STATE AF = 0.; // Flux term A = mu/2(gradUi + gradUiT)(gradUj + gradUjT)
        TPZFNMatrix<9, STATE> DUnj(3,3,0.);
        
        for(int row=0; row<3; row++){
            for(int col=0; col<3; col++){
                DUnj(row, col) = 0.5*(gradUn(row, col) + gradUn(col, row));
            }
        }
        
        AF = TensorInnerProduct(DUi, DUnj);
        ef(i) += 2.*fviscosity*weight*(-AF);
        
        STATE BF = 0.; // Mixed Term B = pdiv(u)
        BF = -pn*divUi;
        ef(i) += weight*(-BF);
        
        for(int j=0; j<nShapeV; j++){
            for(int row=0; row<3; row++){
                phiVj(row,0) = NormalVec(row,j);
            }
            
            TPZFNMatrix<9, STATE> gradVj(3,3,0.);
            TPZFNMatrix<9, STATE> DUj(3,3,0.);
            
            for(int row=0; row<3; row++){
                for(int col=0; col<3; col++){
                    gradVj(row, col) = GradNormVec[j](row, col);
                    gradVj(col, row) = gradVj(row, col);
                    
                    DUj(row, col) = 0.5*(gradVj(row, col) + gradVj(col, row));
                }
            }
            
            STATE Aterm = TensorInnerProduct(DUi, DUj);
            ek(i,j) += 2.0*weight*fviscosity*Aterm; // A - Bilinear gradU*gradV
        }
        
        // Pressure Depending Term
        for(int j=0; j<nShapeP; j++){
            STATE Bterm = -1*weight*phiP(j,0)*divUi;
            
            ek(i, nShapeV+j) += Bterm;
            ek(nShapeV+j, i) += Bterm;
        }
    }
    
    // pablo said to verify this
    for(int i=0; i<nShapeP; i++){
        STATE BtermF = -phiP(i,0)*divsol;
        
        ef(i+nShapeV) += weight*(-BtermF);
    }
    
}

void TPZStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    if(datavec.size() != 2){
        std::cout << "ERROR: DATAVEC SIZE IS DIFFERENT THAN 2\n\n";
        DebugStop();
    }
    
    if(datavec[fVindex].fVecShapeIndex.size()==0) DebugStop();
    
    // Setting the phis
    // for V
    TPZFMatrix<REAL>& phiV = datavec[fVindex].phi;
    TPZFMatrix<REAL>& dphiV = datavec[fVindex].dphix;
    
    TPZFMatrix<REAL>& phiP = datavec[fPindex].phi;
    
    // getting the linear combination or FE approximations
    TPZManVector<STATE> vh = datavec[fVindex].sol[0];
    TPZManVector<STATE> ph = datavec[fPindex].sol[0];
//    TPZFNMatrix<220, REAL> dphiVx(fdimension, dphiV.Cols());
    
    int64_t nShapeV = datavec[fVindex].fDeformedDirections.Cols();
    int64_t nShapeP = phiP.Rows();
    
    int64_t normVecRows = datavec[fVindex].fDeformedDirections.Rows();
    int64_t normVecCols = datavec[fVindex].fDeformedDirections.Cols();
    
    TPZFNMatrix<3, REAL> normVec(normVecRows, normVecCols, 0.);
    normVec = datavec[fVindex].fDeformedDirections;
    
    if(fSpace==1){
        nShapeV /= 2.;
    }
    
    int64_t gy = vh.size();
    
    TPZFNMatrix<9, STATE> phiVi(fdimension, 1, 0.);
    TPZFNMatrix<9, STATE> phiVni(1, 1, 0.);
    
    TPZFNMatrix<9, STATE> phiVj(fdimension, 1, 0.);
    TPZFNMatrix<9, STATE> phiVnj(1, 1, 0.);
    
    TPZFNMatrix<9, STATE> phiPi(fdimension, 1);
    TPZFNMatrix<9, STATE> phiPj(fdimension, 1);
    
    TPZFNMatrix<3,STATE> v1 = bc.Val1();
    TPZManVector<STATE, 3> v2 = bc.Val2();
    STATE pD = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: // Normal Velocity
        {
            if(fSpace==1){
                for(int i=0; i<nShapeV; i++){ // Hdiv adaptation
                    TPZManVector<REAL> n = datavec[fVindex].normal;
                    
                    REAL vhn = vh[0];
                    REAL vn = n[0]*v2[0] + n[1]*v2[1];
                    
                    ef(i,0) += -weight*fBigNumber*(vhn-vn)*phiV(i,0);
                    
                    for(int j=0; j<nShapeV; j++){
                        ek(i,j) += weight*fBigNumber*phiV(j,0)*phiV(i,0);
                    }
                }
                
            } else {
                
                for(int i=0; i<nShapeV; i++){
                    for(int row=0; row<fdimension; row++){
                        phiVi(row,0) = normVec(row, i);
                    }
                    
                    STATE facteF = 0.;
                    for(int is=0; is<gy; is++){
                        facteF += -(vh[is] - v2[is])*phiVi(is,0);
                    }
                    
                    ef(i, 0) += weight*fBigNumber*facteF;
                    
                    for(int j=0; j<nShapeV; j++){
                        for(int row=0; row<fdimension; row++){
                            phiVj(row,0) = normVec(row, j);
                        }
                        
                        STATE facteK = 0.;
                        for(int is=0; is<gy; is++){
                            facteK += phiVj(is,0)*phiVi(is,0);
                        }
                        
                        ek(i,j) += weight*fBigNumber*facteK;
                    }
                }
            }
        }
            break;
            
        case 1: // Tangential Velocity
        {
            DebugStop();
            
            
        }
            break;
            
        case 2: // Normal Stress
        {
            TPZManVector<REAL> n = datavec[fVindex].normal;
            TPZFNMatrix<9, STATE> pn(fdimension, 1);
            
            for(int i=0; i<nShapeV; i++){
                for(int row=0; row<fdimension; row++){
                    phiVi(row,0) = normVec(row, i);
                    pn(row, 0) = n[row]*v1(0,0);
                }
                
                STATE factef = 0.;
                for(int is=0; is<gy; is++){
                    factef += pn(is,0)*phiVi(is,0);
                }
                
                ef(i,0) += weight*factef;
            }
        }
            break;
            
        case 3: // Tangential Stress
        {
            DebugStop();        }
            break;
            
        default:
        {
            std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
            DebugStop();
        }
            break;
    }
    
    
}

int TPZStokesMaterial::VariableIndex(const std::string& name) const {
    
    if(!strcmp("P", name.c_str())) return 0;
    if(!strcmp("Pressure", name.c_str())) return 0;
    if(!strcmp("State", name.c_str())) return 0;
    if(!strcmp("V", name.c_str())) return 1;
    if(!strcmp("f", name.c_str())) return 2;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZStokesMaterial::NSolutionVariables(int var) const{
    
    int aux;
    switch (var) {
        case 0: // pressure  [scalar]
            aux = 1;
            break;
            
        case 1: // velocity [vector]
        case 2: // external force [vector]
            aux = 3;
            break;
            
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
    }
    return aux;
}

void TPZStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    TPZManVector<STATE, 3> v_h = datavec[fVindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[fPindex].sol[0];
    
    TPZFNMatrix<9, STATE> gradU(3,1);
    
    TPZFMatrix<STATE>& dsol = datavec[fVindex].dsol[0];
    TPZFNMatrix<9,STATE> dsolxy(3,3), dsolxyp(3,1);
    
    dsolxy = dsol;
    
    Solout.Resize(NSolutionVariables(var));
    
    switch(var) {
            
        case 0:{
            Solout[0] = p_h[0];
            break;
        }
            
        case 1:{
            Solout[0] = v_h[0]; //Vx
            Solout[1] = v_h[1]; //Vy
            Solout[2] = v_h[2]; //Vz
            break;
        }
            
        case 2:{
            TPZVec<STATE> f(3,0.);
            
            if(this->HasForcingFunction()){
                this->ForcingFunction()(datavec[fVindex].x, f);
            }
            
            Solout[0] = f[0];
            Solout[1] = f[1];
            Solout[2] = f[2];
            break;
        }
            
        default:{
            std::cout << "\n\nVar index not implemented\n\n";
            DebugStop();
        }
    }
}

void TPZStokesMaterial::ContributeInterface(const TPZMaterialDataT<STATE>& data, const std::map<int, TPZMaterialDataT<STATE>>& datavecleft, const std::map<int, TPZMaterialDataT<STATE>>& datavecright, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef){
    
    int nRefLeft = datavecleft.size();
    int nRefRight = datavecright.size();
    
    if(datavecleft.find(fVindex) == datavecleft.end()) DebugStop();
    if(datavecright.find(fPindex) == datavecright.end()) DebugStop();
    
    const TPZMaterialDataT<STATE>& vDataLeft = datavecleft.find(fVindex)->second;
    const TPZMaterialDataT<STATE>& pDataRight = datavecright.find(fPindex)->second;
    
    if(vDataLeft.fVecShapeIndex.size() == 0) DebugStop();
    
    // Setting the phis
    // V -> Left
    const TPZFNMatrix<9, REAL>& tan = pDataRight.axes;
    TPZManVector<STATE, 3> un = vDataLeft.sol[0];
    STATE sLambdan = pDataRight.sol[0][0];
    TPZManVector<STATE, 3> Lambdan = pDataRight.sol[0];
        
    int nShapeV = vDataLeft.fVecShapeIndex.NElements();
    int nShapeLambda = pDataRight.phi.Rows();
    int nStateVariablesL = fdimension-1;
    
    int normVecRows = vDataLeft.fDeformedDirections.Rows();
    int normVecCols = vDataLeft.fDeformedDirections.Cols();
    TPZFNMatrix<3, REAL> NormalVec(normVecRows, normVecCols, 0.);
    
    if(vDataLeft.fNeedsDeformedDirectionsFad){
        for(int row=0; row<normVecRows; row++){
            for(int col=0; col<normVecCols; col++){
                NormalVec(row, col) = vDataLeft.fDeformedDirectionsFad(row,col).val();
            }
        }
        
    } else {
        NormalVec = vDataLeft.fDeformedDirections;
    }
    
    for(int i=0; i<nShapeV; i++){
        TPZFNMatrix<9, STATE> phiVi(3,1,0.);
        TPZFNMatrix<9, STATE> lambdan(3,1,0.);
        
        for(int row=0; row<3; row++){
            phiVi(row,0) = NormalVec(row, i);
        }
        
        STATE LambdaDotPhiV = 0.;
        for(int f=0; f<nStateVariablesL; f++){
            for(int row=0; row<3; row++){
                
            }
        }
    }
    
}

void TPZStokesMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                      const std::map<int, TPZMaterialDataT<STATE>> &datavec,
                      REAL weight,
                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                       TPZBndCondT<STATE> &bc){
    DebugStop();
}

void TPZStokesMaterial::FillDataRequirementsInterface(TPZMaterialDataT<STATE>& data, std::map<int, TPZMaterialDataT<STATE>>& datavec_left, std::map<int, TPZMaterialDataT<STATE>>& datavec_right){
    int nref_left = datavec_left.size();
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsSol = true;
    datavec_right[1].fNeedsNormal = true;
    datavec_right[1].fNeedsSol = true;
    datavec_left[0].fNeedsDeformedDirectionsFad = true;
}
