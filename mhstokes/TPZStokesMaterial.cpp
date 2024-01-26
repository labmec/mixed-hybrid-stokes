#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>

#include "TPZStokesMaterial.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesmaterial");
#endif

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

    int64_t dimension = Dimension(); // problems dimension

    int64_t nShapeV = datavec[EVindex].fVecShapeIndex.NElements(); // number of velocity Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiV(dimension, nShapeV, 0.0);
    TPZFNMatrix<20, REAL>& divPhiV = datavec[EVindex].divphi;

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

    const int nterms = dimension * (dimension + 1) / 2;
    TPZFNMatrix<150, REAL> StrainRate(nterms, nShapeV, 0.0); //Using voight notation
    
    if (datavec[EVindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {
            int cont = dimension-1;
            for (int64_t i = 0; i < dimension; i++)
            {
                PhiV(i, j) = datavec[EVindex].fDeformedDirectionsFad(i, j).val();
                StrainRate(i,j) = datavec[EVindex].fDeformedDirectionsFad(i, j).fastAccessDx(i);
                for (int64_t k = i+1; k < dimension && k != i; k++)
                {
                    StrainRate(++cont,j) = 0.5*(datavec[EVindex].fDeformedDirectionsFad(i, j).fastAccessDx(k) + datavec[EVindex].fDeformedDirectionsFad(k, j).fastAccessDx(i));
                }
            }
        }
    }

    TPZFNMatrix<9, REAL> VoightCorrection(nterms, nterms, 0.0); //This is to multiply by 2 the off diagonal part of strain rate tensor to account for its symmetry
    for (int64_t i = 0; i < nterms; i++)
    {
        VoightCorrection(i,i) = i<dimension ? 1.0 : 2.0;
    }
    
    TPZFNMatrix<3,REAL> SourceTerm(dimension, 1.0, 0.0);
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EVindex].x, sourceAux);
        for (int64_t i = 0; i < dimension; i++)
        {
            SourceTerm(i,0) = sourceAux[i] * fviscosity;
        }
    }

    //Body Forces contribution
    ef.AddContribution(0, 0, PhiV, true, SourceTerm, false, weight);

    //Flux Matrix A contribution
    TPZFNMatrix<150, REAL> matrixC;
    VoightCorrection.Multiply(StrainRate, matrixC);
   
    REAL factor = 2.0 * fviscosity * weight;
    ek.AddContribution(0, 0, StrainRate, true, matrixC, false, factor);
    
    factor = fviscosity*weight;
    ek.AddContribution(0, 0, divPhiV, false, divPhiV, true, factor);

    //Divergence Matrix B contribution
    factor = -1.0 * weight;
    ek.AddContribution(0, nShapeV, divPhiV, false, PhiP, true, factor);

    //Divergence Matrix BT contribution
    ek.AddContribution(nShapeV, 0, PhiP, false, divPhiV, true, factor);
    
    if(datavec.size()>2)
    {
        TPZFMatrix<REAL>& phiG = datavec[EGindex].phi;
        TPZFMatrix<REAL>& phipM = datavec[EPMindex].phi;
        
        // Pressure and distributed flux
        for(int j=0; j<nShapeP; j++)
        {
            ek(nShapeV+nShapeP, nShapeV+j) += PhiP(j,0)*phiG(0,0)*weight;
            ek(nShapeV+j, nShapeV+nShapeP) += PhiP(j,0)*phiG(0,0)*weight;
        }
        
        // Injection and average-pressure
        ek(nShapeV+nShapeP+1, nShapeV+nShapeP) += phiG(0,0)*phipM(0,0)*weight;
        ek(nShapeV+nShapeP, nShapeV+nShapeP+1) += phiG(0,0)*phipM(0,0)*weight;
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

void TPZStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    TPZFNMatrix<150, REAL> PhiV = datavec[EVindex].phi;
    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    
    int64_t nShapeV = PhiV.Rows();
    int64_t nShapeP = PhiP.Rows();
            
    TPZFNMatrix<20,STATE> val1(3,3,0.); //grad
    TPZManVector<STATE, 4> val2(3, 0.); //value
    REAL pressure = 0.0;
    
    if (bc.HasForcingFunctionBC())
    {
        TPZVec<STATE> vVal;
        TPZFMatrix<STATE> gradval;
        bc.ForcingFunctionBC()(datavec[EVindex].x,val2,val1);
        pressure = val2[3];
    }
    else
    {
        val1 = bc.Val1();
        val2 = bc.Val2();
    }
    
    switch (bc.Type()) {
        case 0: // Normal Velocity
        {
            REAL v_n = val2[0]; // if bc was set in .json file, the normal value was already pescribed
            if (bc.HasForcingFunctionBC()) // if the bc is set through an analytic solution, we need to compute its normal component
            {
                v_n = 0.0;
                for (int i=0; i<fdimension; i++)
                {
                    v_n += val2[i]*datavec[EVindex].normal[i];
                }
            }
            
            REAL factor = fBigNumber * weight;
            
            for(int64_t j = 0; j < nShapeV; j++)
            {
                ef(j) += v_n * PhiV(j,0) * factor;
                
                for(int64_t i = 0; i < nShapeV; i++)
                {
                    ek(i, j) += PhiV(i, 0) * PhiV(j, 0) * factor;
                }
            }
            break;
        }
            
        case 1: // Tangential Velocity
        {
            TPZManVector<REAL, 3> v_t = {0., 0., 0.}; // for tangential bc, a vector is prescribed, so we take the inner product with the local tangential axe
            for(int i = 0; i < fdimension-1; i++){
                for(int j = 0; j < fdimension; j++){
                    v_t[i] += val2[j] * datavec[EVindex].axes(i,j);
                }
            }
            
            REAL factor = fBigNumber * weight;
            
            for (int64_t j = 0; j < nShapeV; j++)
            {
                for (int64_t k = 0; k < fdimension-1; k++)
                {
                    int64_t index1 = (fdimension-1) * j + k;
                    ef(index1) += -v_t[k] * PhiV(j, 0) * factor;
                    
                    for (int64_t i = 0; i < nShapeV; i++)
                    {
                        for (int64_t l = 0; l < fdimension-1; l++)
                        {
                            int64_t index2 = (fdimension-1) * i + l;
                            if (k != l) continue;
                            ek(index1, index2) += PhiV(i, 0)*PhiV(j, 0)*factor;
                        }
                    }
                }
            }
            break;
        }
            
        case 2: // Normal Stress
        {
            REAL sigma_nn = val2[0]; // if bc was set in .json file, the normal value was already prescribed
            if(bc.HasForcingFunctionBC()){ // if the bc is set through an analytic solution, we need to compute its normal component from the velocity gradient
                const int n = fdimension * (fdimension+1)/2; //size of the vector in Voight notation
                
                TPZFNMatrix<6, REAL> sigmaVoight(n, 1, 0.);
                StressTensor(val1, sigmaVoight, pressure);
                
                TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
                int cont = fdimension-1;
                
                for(int i = 0; i < fdimension; i++){
                    sigma(i, i) = sigmaVoight(i, 0);
                    
                    for(int j = i + 1; j < fdimension; j++){
                        sigma(i, j) = sigmaVoight(++cont,0);
                        sigma(j, i) = sigmaVoight(cont, 0);
                    }
                }
                
                TPZFNMatrix<3, REAL> sigma_n(fdimension, 1, 0.0);
                
                for(int i = 0; i < fdimension; i++){
                    for(int j = 0; j < fdimension; j++){
                        sigma_n(i, 0) += sigma(i, j)*datavec[EVindex].normal[j];
                    }
                }
                
                sigma_nn = 0.0;
                for(int i = 0; i < fdimension; i++){
                    sigma_nn += sigma_n[i]*datavec[EVindex].normal[i];
                }
            }
            
            for(int64_t i = 0; i < nShapeV; i++){
                ef(i) += sigma_nn * PhiV(i, 0) * weight;
            }
        }
            break;
            
        case 3: // Tangential Stress
        {
            TPZManVector<REAL,3> sigma_nt(fdimension - 1, 0.);
            if(bc.HasForcingFunctionBC()){
                const int n = fdimension * (fdimension - 1) / 2;
                
                TPZFNMatrix<6, REAL> sigmaVoight(n, 1, 0.0);
                StressTensor(val1, sigmaVoight, pressure);
                
                TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
                int cont = fdimension-1;
                
                for(int  i=0; i<fdimension; i++){
                    sigma(i,i) = sigmaVoight(i,0);
                    
                    for(int j=0; j<fdimension; j++){
                        sigma(i,j) = sigmaVoight(++cont,0);
                        sigma(j,i) = sigmaVoight(cont,0);
                    }
                }
                
                TPZManVector<REAL, 3> sigma_n(fdimension,0.0);
                
                for(int i=0; i<fdimension; i++){
                    for(int j=0; j<fdimension; j++){
                        sigma_n[i] += sigma(i,j)*datavec[EVindex].normal[j];
                    }
                }
                
                for(int i=0; i<fdimension; i++){
                    for(int j=0; j<fdimension; j++){
                        sigma_nt[i] += sigma_n[i]*datavec[EVindex].axes(i,j);
                    }
                }
            }
            else
            {
                for (int i = 0; i < fdimension-1; i++)
                    for (int j = 0; j < fdimension; j++)
                        sigma_nt[i] += val2[j]*datavec[EVindex].axes(i,j);
            }
            
            for (int64_t i = 0; i < nShapeV; i++)
            {
                for (int j = 0; j < fdimension-1; j++)
                {
                    int64_t index = (fdimension-1) * i + j;
                    ef(index) += -PhiV(i, 0) * sigma_nt[j] * weight;
                }
            }
            
            break;
        }
            
        default:
        {
            std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
            DebugStop();
            break;
        }
    }
}

int TPZStokesMaterial::VariableIndex(const std::string& name) const {
    
    if(!strcmp("Pressure", name.c_str())) return 0;
    if(!strcmp("Velocity", name.c_str())) return 1;
    if(!strcmp("Force", name.c_str())) return 2;
    if(!strcmp("Stress", name.c_str())) return 3;
    
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
        case 3: // stress tensor
            aux = 9;
            break;
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
    }
    return aux;
}

void TPZStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    TPZManVector<STATE, 3> v_h = datavec[EVindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];
    
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
                this->ForcingFunction()(datavec[EVindex].x, f);
            }
            
            Solout[0] = f[0];
            Solout[1] = f[1];
            Solout[2] = f[2];
            break;
        }
        case 3: {
            
                TPZFNMatrix<10,STATE> gradUn = datavec[EVindex].dsol[0];
                TPZFNMatrix<9,STATE> DUn_j(3,3,0.), sigma(3,3,0.);
            
                for (int e=0; e<Dimension(); e++) {
                    for (int f=0; f<Dimension(); f++) {
                         DUn_j(e,f)= 0.5 * (gradUn(e,f) + gradUn(f,e));
                        sigma(e,f) = 2.*fviscosity*DUn_j(e,f);
                    }
                    sigma(e,e) -= p_h[0];
                }
            
                for (int e=0; e<3; e++) {
                    for (int f=0; f<3; f++) {
                        Solout[e*3+f] = sigma(e,f);
                    }
                }
            }
            
            break;
            
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

void TPZStokesMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    datavec[EVindex].fNeedsSol = false;
    datavec[EPindex].fNeedsSol = false;
    datavec[EVindex].fNeedsNormal = true;
    datavec[EPindex].fNeedsNormal = true;
}

void TPZStokesMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors){
    
    if(!HasExactSol()) DebugStop();
    
    TPZManVector<STATE, 4> sol_exact(3);
    TPZFNMatrix<9> dsol_exact(3,3);
    
    fExactSol(data[EVindex].x, sol_exact, dsol_exact);
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.);
    
    TPZManVector<STATE> Velocity(3,0.);
    TPZManVector<STATE> Pressure(3,0.);
    
    this->Solution(data, VariableIndex("Velocity"), Velocity);
    this->Solution(data, VariableIndex("Pressure"), Pressure);
    
    TPZFMatrix<STATE>& dsolv = data[EVindex].dsol[0];
    
    dsolv.Resize(3, 3);
    
    STATE diffv, diffp, diffdiv;
    
    errors[0] = 0.;
    for(int i=0; i<3; i++){
        diffv = Velocity[i] - sol_exact[i];
        errors[0] += diffv*diffv;
    }
    
    STATE div_exact=0., Div=0.;
    for(int i=0; i<3; i++){
        div_exact += dsol_exact[i];
        Div += dsolv(i, i);
    }
    
    diffdiv = Div - div_exact;
    errors[1] = diffdiv*diffdiv;
    
    diffp = Pressure[0] - sol_exact[3];
    errors[2] = diffp*diffp;
}

void TPZStokesMaterial::StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma, REAL pressure){
    
    const int n = fdimension*(fdimension+1)/2;
    TPZFNMatrix<6, REAL> strain(n,1,0.0);
    TPZFNMatrix<36, REAL> D(n,n,0.0);
    
    StrainTensor(gradU, strain);
    ViscosityTensor(D);
    
    D.Multiply(strain, sigma);
    
    for(int i=0; i<fdimension; i++){
        sigma(i,0) -= pressure;
    }
}

void TPZStokesMaterial::StrainTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6, REAL>& epsilon){
    
    const int n = fdimension*(fdimension+1)/2;
    int cont = fdimension-1;
    
    for(int i=0; i<fdimension; i++){
        
        epsilon(i,0) = gradU(i,i);
        
        for(int j=i+1; j<fdimension; j++){
            epsilon(++cont, 0) = 0.5*(gradU(i,j)+gradU(j,i));
        }
    }
}

void TPZStokesMaterial::ViscosityTensor(TPZFNMatrix<36, REAL>& D){
    int n=fdimension*(fdimension+1)/2;
    
    for(int i=0; i<n; i++){
        D(i,i) = 2*fviscosity;
    }
}

//void TPZStokesMaterial::ToVoigt(const TPZFNMatrix<9, STATE> &Sigma, TPZFNMatrix<6, STATE>& Svoigt) const {
//    Svoigt(Exx,0) = Sigma(0,0);
//    Svoigt(Exy,0) = Sigma(0,1);
//    Svoigt(Eyy,0) = Sigma(1,1);
//
//    if(fdimension>2){
//        Svoigt(Exz,0) = Sigma(0,2);
//        Svoigt(Eyz,0) = Sigma(1,2);
//        Svoigt(Ezz,0) = Sigma(2,2);
//    }
//}
//
//void TPZStokesMaterial::FromVoigt(const TPZFNMatrix<6,STATE> &Svoigt, TPZFNMatrix<9, STATE>& Sigma) const {
//    Sigma(0,0) = Svoigt(Exx, 0);
//    Sigma(0,1) = Svoigt(Exy, 0);
//    Sigma(1,1) = Svoigt(Eyy, 0);
//
//    if(fdimension > 2){
//        Sigma(0,2) = Svoigt(Exz,0);
//        Sigma(1,2) = Svoigt(Eyz,0);
//        Sigma(2,2) = Svoigt(Ezz,0);
//    }
//}
