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
    if (datavec[fVindex].fVecShapeIndex.size() == 0) {
        DebugStop();
    }
    
    NormalVec = datavec[fVindex].fDeformedDirections;
    
    // Setting the phis (is this the correct variable? or should it be fDeformedDirections?)
    
    
    
    // for P
    TPZFMatrix<REAL>& phiP = datavec[fPindex].phi;
    
    int64_t nShapeV = datavec[fVindex].fVecShapeIndex.NElements();
    int64_t nShapeP = phiP.Rows();
    
    // creating a vector of gradients of the H(div) vectors (why is it that way?)
    TPZManVector<TPZFNMatrix<9, REAL>, 18> GradNormVec(normVecCols);
    
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
                for(int j=0; j < this->Dimension(); j++){
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
            }
        }
        
        for(int row=0; row<3; row++){
            for(int col=0; col<3; col++){
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
        
        for(int j=0; j<nShapeV; j++){
            for(int row=0; row<3; row++){
                phiVj(row,0) = NormalVec(row,j);
            }
            
            TPZFNMatrix<9, STATE> gradVj(3,3,0.);
            TPZFNMatrix<9, STATE> DUj(3,3,0.);
            
            for(int row=0; row<3; row++){
                for(int col=0; col<3; col++){
                    gradVj(row, col) = GradNormVec[j](row, col);
                }
            }
            
            for(int row=0; row<3; row++){
                for(int col=0; col<3; col++){
                    // symetric gradient of the test function
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
    
    ek.Print("ek",std::cout,EMathematicaInput);
    ef.Print("ef",std::cout,EMathematicaInput);
}

void TPZStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    if(datavec.size() != 2){
        std::cout << "ERROR: DATAVEC SIZE IS DIFFERENT THAN 2\n\n";
        DebugStop();
    }
    
//    if(datavec[fVindex].fVecShapeIndex.size()==0) DebugStop();
    
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
//            DebugStop();
            
        }
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

void TPZStokesMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsHSize = true;
        datavec[idata].fNeedsNormal = true;
    }
    datavec[0].fNeedsDeformedDirectionsFad = true;
}
