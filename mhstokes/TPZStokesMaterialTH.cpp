#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>
#include <fstream>
#include <ostream>
#include <iostream>

#include "TPZStokesMaterialTH.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesmaterial");
#endif

TPZStokesMaterialTH::TPZStokesMaterialTH(): TBase() {
    
}

TPZStokesMaterialTH::TPZStokesMaterialTH(int matID, int dimension, double viscosity) : TBase(matID), fdimension(dimension), fviscosity(viscosity){
    
}

TPZStokesMaterialTH::~TPZStokesMaterialTH(){
    
}

template <typename TVar>
TVar TPZStokesMaterialTH::TensorInnerProduct(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    if(S.Rows() != S.Cols() || T.Rows() != T.Cols() || S.Rows() != T.Rows()) DebugStop();
    
    TVar Val = 0;
    
    for(int i=0; i<S.Cols(); i++){
        for(int j=0; j<S.Rows(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
}

void TPZStokesMaterialTH::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef)
{
    STATE mu = fviscosity;
    
    TPZFMatrix<REAL> &PhiV = datavec[EVindex].phi;
    int64_t nShapeV = PhiV.Rows();
    
    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows();
    
    TPZFNMatrix<60, REAL> dphi = datavec[EVindex].dphix;
    auto axes = datavec[EVindex].axes;
    
    TPZFNMatrix<3, REAL> dphiV(fdimension, nShapeV, 0.0);
    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphiV, axes);
    
    TPZFNMatrix<3, REAL> SourceTerm(fdimension, 1, 0.0);
    TPZVec<REAL> sourceAux(3);
    
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EVindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
            SourceTerm(i, 0) = sourceAux[i];
    }
    
    const int dim = Dimension();
    int n = dim * (dim + 1) / 2;
    
    // divergence of v
    TPZFNMatrix<150, STATE> divPhiV(dim * nShapeV, 1, 0.0);
    for (int j = 0; j < nShapeV; j++)
        for (int i = 0; i < dim; i++)
            divPhiV(j * dim + i, 0) = dphiV(i, j);
    
    // strain rate tensor
    TPZFNMatrix<150, STATE> strainRate(n, nShapeV * dim, 0.0);

    for (int j = 0; j < nShapeV; j++)
    {
        int cont = dim;
        for (int i = 0; i < dim; i++)
        {
            strainRate(i, j * dim + i) = dphiV(j, i);
            
            for (int k = i + 1; k < dim; k++)
            {
                strainRate(cont, dim * j + i) = dphiV(j, k);
                strainRate(cont, dim * j + k) = dphiV(j, i);
                cont++;
            }
        }
    }
    
    TPZFMatrix<REAL> phiV_force(dim * nShapeV, dim, 0.0);
    for (int j = 0; j < nShapeV; j++)
    {
        for (int i = 0; i < dim; i++)
        {
            phiV_force(dim * j + i, i) = PhiV(j);
        }
    }
    
    // this is to multiply by 2 the off diagonal part of strain tensor
    TPZFNMatrix<9, REAL> VoigtCorrection(n, n, 0.0);
    for (int64_t i = 0; i < n; i++)
        VoigtCorrection(i, i) = (i < dim) ? 1.0 : 0.5;
    
    // body forces contribution
    ef.AddContribution(0, 0, phiV_force, false, SourceTerm, false, weight);
    
    // flux - Matrix A contribution
    TPZFNMatrix<150, REAL> matrixC;
    VoigtCorrection.Multiply(strainRate, matrixC);
    
    STATE factor = 2.0 * mu * weight;
    ek.AddContribution(0, 0, strainRate, true, matrixC, false, factor);
    
    // Divergence - Matrix B contribution
    factor = - 1.0 * weight;
    ek.AddContribution(0, nShapeV, divPhiV, false, PhiP, true, factor);
    
    // Divergence - Matrix BT contribution
    ek.AddContribution(nShapeV, 0, PhiP, false, divPhiV, true, factor);
 
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

void TPZStokesMaterialTH::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    TPZFNMatrix<150, REAL> phiV = datavec[EVindex].phi;
    TPZFMatrix<REAL> &phiP = datavec[EPindex].phi;
    
    int64_t nShapeV = phiV.Rows();
    int64_t nShapeP = phiP.Rows();
    
    TPZFNMatrix<20, STATE> val1(3, 3, 0.0);
    TPZManVector<STATE, 3> val2(3, 0.0);
    const auto &BIGNUMBER = fBigNumber;
    
    if (bc.HasForcingFunctionBC())
    {
        bc.ForcingFunctionBC()(datavec[EVindex].x, val2, val1);
    }
    else
    {
        val1 = bc.Val1();
        val2 = bc.Val2();
    }
    
    int dim = fdimension;
    
    switch (bc.Type())
    {
        case 0: // Dirichlet condition
        {
            for (int j = 0; j < nShapeV; j++)
            {
                for (int i = 0; i < dim; i++)
                {
                    ef(dim * j + i, 0) += phiV(j, 0) * val2[i] * BIGNUMBER * weight;
                    
                    for (int k = 0; k < nShapeV; k++)
                    {
                        ek(dim * j + i, dim * k + i) += phiV(j, 0) * phiV(k, 0) * BIGNUMBER * weight;
                    }
                }
            }
        }
            break;
        
        case 1: // Neumann condition
        {
            for (int j = 0; j < nShapeV; j++)
            {
                for (int i = 0; i < dim; i++)
                {
                    ef(dim * j + i, 0) += val2[i] * phiV(j, 0) * weight;
                }
            }
        }
            break;
            
        case 2: // Mixed Condition
        {
            DebugStop(); // please implement me
        }
            break;
            
        default:
        {
            std::cout << "ERROR: BOUNDARY CONIDITON NOT IMPLEMENTED" << std::endl;
            DebugStop();
        }
            break;
    }
}

int TPZStokesMaterialTH::VariableIndex(const std::string& name) const {
    
    // Pressure
    if (!strcmp("Pressure", name.c_str()))
        return EPressure;
    if (!strcmp("ExactPressure", name.c_str()))
        return EExactPressure;
    if (!strcmp("ErrorPressure", name.c_str()))
        return EErrorPressure;
    
    // Velocity
    if (!strcmp("Velocity", name.c_str()))
        return EVelocity;
    if (!strcmp("ExactVelocity", name.c_str()))
        return EExactVelocity;
    if (!strcmp("ErrorVelocity", name.c_str()))
        return EErrorVelocity;
    
    // Source term
    if (!strcmp("SourceTerm", name.c_str()))
        return ESourceTerm;
    
    // Stress
    if (!strcmp("Stress", name.c_str()))
        return EStress;
    if (!strcmp("ExactStress", name.c_str()))
        return EExactStress;
    if (!strcmp("ErrorStress", name.c_str()))
        return EErrorStress;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZStokesMaterialTH::NSolutionVariables(int var) const{
    
    int aux;
    
    switch (var) {
        case EPressure:
        case EExactPressure:
        case EErrorPressure:
            aux = 1;
            break;
            
        case EVelocity:
        case EExactVelocity:
        case EErrorVelocity:
        case ESourceTerm:
            aux = 3;
            break;
        
        case EStress:
        case EExactStress:
        case EErrorStress:
            aux = 9;
            break;
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
    }
    return aux;
}

void TPZStokesMaterialTH::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    const int n = fdimension * (fdimension + 1) / 2;
    
    TPZManVector<STATE, 3> u_h = datavec[EVindex].sol[0];
    TPZFNMatrix<10,STATE> gradU_h = datavec[EVindex].dsol[0];
    
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];
    
    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9, STATE> gradsol_exact(3, 3);
    STATE p_exact = 0.;
    
    if (this->HasExactSol())
    {
        fExactSol(datavec[EVindex].x, sol_exact, gradsol_exact);
        p_exact = sol_exact[3];
    }
    
    Solout.Resize(NSolutionVariables(var));
    
    switch(var) {
            
        case EPressure: // pressure
        {
            Solout[0] = p_h[0];
        }
            break;
            
        case EExactPressure: // exact pressure
        {
            Solout[0] = p_exact;
        }
            break;
        
        case EErrorPressure: // pressure element error
        {
            Solout[0] = abs(p_exact - p_h[0]);
        }
            break;
            
        case EVelocity: // velocity
        {
            Solout[0] = u_h[0]; // vx
            Solout[1] = u_h[1]; // vy
            Solout[2] = u_h[2]; // vz
        }
            break;
        
        case EExactVelocity: // exact velocity
        {
            Solout[0] = sol_exact[0];
            Solout[1] = sol_exact[1];
            Solout[2] = sol_exact[2];
        }
            break;

        case EErrorVelocity: // velocity element error
        {
            Solout[0] = abs(sol_exact[0] - u_h[0]); // vx
            Solout[1] = abs(sol_exact[1] - u_h[1]); // vy
            Solout[2] = abs(sol_exact[2] - u_h[2]); // vz
        }
            break;
            
        case ESourceTerm: // source term
        {
            TPZManVector<STATE, 3> f(3, 0.);
            
            if (!this->HasForcingFunction())
                this->ForcingFunction()(datavec[EVindex].x, f);
            
            Solout[0] = f[0]; // fx
            Solout[1] = f[1]; // fy
            Solout[2] = f[2]; // fz
        }
            break;
            
        case EStress: // stress
        {
            TPZFNMatrix<6, STATE> sigmaVoigt(n, 1, 0.);
            TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
            
            StressTensor(gradU_h, sigmaVoigt, p_h[0]);
            FromVoigt(sigmaVoigt, sigma);
            
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Solout[i * 3 + j] = sigma(i, j);
        }
            break;
            
        case EExactStress: // exact stress
        {
            TPZFNMatrix<6, STATE> sigmaVoigt_exact(n, 1, 0.0);
            TPZFNMatrix<9, STATE> sigma_exact(3, 3, 0.0);

            StressTensor(gradsol_exact, sigmaVoigt_exact, p_exact);
            FromVoigt(sigmaVoigt_exact, sigma_exact);
            
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Solout[i * 3 + j] = sigma_exact(i, j);
        }
            break;
            
        case EErrorStress: // stress element error
        {   
            TPZFNMatrix<6, STATE> sigmaVoigt(n, 1, 0.0), sigmaVoigt_exact(n, 1, 0.0);
            TPZFNMatrix<9, STATE> sigma_h(3, 3, 0.0), sigma_exact(3, 3, 0.0);
            
            StressTensor(gradU_h, sigmaVoigt, p_h[0]);
            StressTensor(gradsol_exact, sigmaVoigt_exact, p_exact);
            
            FromVoigt(sigmaVoigt, sigma_h);
            FromVoigt(sigmaVoigt_exact, sigma_exact);
            
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Solout[i * 3 + j] = abs(sigma_exact(i, j) - sigma_h(i, j));
        }
            break;
            
        default:{
            std::cout << "\n\nVar index not implemented\n\n";
            DebugStop();
        }
    }
}

void TPZStokesMaterialTH::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
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

void TPZStokesMaterialTH::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    datavec[EVindex].fNeedsSol = false;
    datavec[EPindex].fNeedsSol = false;
    datavec[EVindex].fNeedsNormal = true;
    datavec[EPindex].fNeedsNormal = true;
}

void TPZStokesMaterialTH::Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors)
{
    // 0: L2 p, 1: L2 p_ex, 2: L2 u, 3: L2 u_ex, 4: L2 divu, 5: L2 divu, 6: L2 sigma, 7: L2 sigma_ex
    
    if (!HasExactSol()) DebugStop();
    
    errors.Resize(NEvalErrors());
    
    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9, STATE> gradsol_exact(3, 3);
    
    // Getting the exact solution for velocity, pressure and velocity gradient
    fExactSol(data[EVindex].x, sol_exact, gradsol_exact);
    STATE p_exact = sol_exact[3];
    
    // Getting the numeric solution for velocity, pressure, and velocity gradient
    TPZManVector<STATE> u_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10, STATE> gradv_h = data[EVindex].dsol[0];
    
    this->Solution(data, VariableIndex("Velocity"), u_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);
    
    STATE diff_v, diff_p, diff_div;
    
    // Pressure error
    diff_p = p_h[0] - p_exact;
    errors[0] = diff_p * diff_p;
    errors[1] = p_exact * p_exact;
    
    // Velocity Error
    errors[2] = 0.0;
    errors[3] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        diff_v = u_h[i] - sol_exact[i];
        errors[2] += diff_v * diff_v;
        errors[3] += sol_exact[i] * sol_exact[i];
    }
    
    // Velocity Divergence Error
    STATE div_exact = 0.0, div_h = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        div_exact += gradsol_exact(i, i);
        div_h += gradv_h(i, i);
    }
    
    diff_div = div_h - div_exact;
    errors[4] = diff_div * diff_div;
    errors[5] = div_exact * div_exact;
    
    // Stress Tensor Error
    const int n = fdimension * (fdimension + 1) / 2;
    TPZFNMatrix<6, REAL> sigma_exact(n, 1), sigma_h(n, 1);
    
    StressTensor(gradsol_exact, sigma_exact, p_exact);
    StressTensor(gradv_h, sigma_h, p_h[0]);
    
    errors[6] = 0.0;
    errors[7] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        const STATE diff_sigma = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += diff_sigma * diff_sigma;
        errors[7] += sigma_exact(i, 0) * sigma_exact(i, 0);
    }
    
    for (int i = fdimension; i < n; i++)
    {
        const STATE diff_sigma = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += 2.0 * diff_sigma * diff_sigma;
        errors[7] += 2.0 * sigma_exact(i, 0) * sigma_exact(i, 0);
    }
    
    // Deviatoric Stress Tensor
    TPZFNMatrix<6, REAL> devSigma_exact(n, 1), devSigma_h(n, 1);
    
    DeviatoricStressTensor(gradsol_exact, devSigma_exact);
    DeviatoricStressTensor(gradv_h, devSigma_h);
    
    errors[8] = 0.0;
    errors[9] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        const STATE diff_devsigma = devSigma_h(i, 0) - devSigma_exact(i, 0);
        errors[8] += diff_devsigma * diff_devsigma;
        errors[9] += devSigma_exact(i, 0) * devSigma_exact(i, 0);
    }
    
    for (int i = fdimension; i < n; i++)
    {
        const STATE diff_devsigma = devSigma_h(i, 0) - devSigma_exact(i, 0);
        errors[8] += 2.0 * diff_devsigma * diff_devsigma;
        errors[9] += 2.0 * devSigma_exact(i, 0) * devSigma_exact(i, 0);
    }
}

void TPZStokesMaterialTH::StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma, REAL pressure){
    
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

void TPZStokesMaterialTH::DeviatoricStressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma){
    
    const int n = fdimension*(fdimension+1)/2;
    TPZFNMatrix<6, REAL> strain(n,1,0.0);
    TPZFNMatrix<36, REAL> D(n,n,0.0);
    
    StrainTensor(gradU, strain);
    ViscosityTensor(D);
    
    D.Multiply(strain, sigma);
}

void TPZStokesMaterialTH::StrainTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6, REAL>& epsilon){
    
    const int n = fdimension*(fdimension+1)/2;
    int cont = fdimension-1;
    
    for(int i=0; i<fdimension; i++){
        
        epsilon(i,0) = gradU(i,i);
        
        for(int j=i+1; j<fdimension; j++){
            epsilon(++cont, 0) = 0.5*(gradU(i,j)+gradU(j,i));
        }
    }
}

void TPZStokesMaterialTH::ViscosityTensor(TPZFNMatrix<36, REAL>& D){
    int n=fdimension*(fdimension+1)/2;
    
    for(int i=0; i<n; i++){
        D(i,i) = 2*fviscosity;
    }
}

void TPZStokesMaterialTH::FromVoigt(const TPZFNMatrix<6,STATE> &sigmaVoigt, TPZFNMatrix<9, STATE>& sigma) const
{
    int cont = fdimension - 1;
    
    for (int i = 0; i < fdimension; i++)
    {
        sigma(i, i) = sigmaVoigt(i, 0);
        
        for (int j = i + 1; j < fdimension; j++)
        {
            sigma(i, j) = sigmaVoigt(++cont, 0);
            sigma(j, i) = sigmaVoigt(cont, 0);
        }
    }
}
