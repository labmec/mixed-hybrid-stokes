#include "TPZAxisymStokesMaterial.h"
#include <pzlog.h>

#ifdef PZ_LOG
static TPZLogger logger("pz.axisymmetricstokesmaterial");
#endif

TPZAxisymStokesMaterial::TPZAxisymStokesMaterial(): TPZStokesMaterial() {}

TPZAxisymStokesMaterial::TPZAxisymStokesMaterial(int matID, int dimension, double viscosity) : TPZStokesMaterial(matID, dimension, viscosity) {}

TPZAxisymStokesMaterial::~TPZAxisymStokesMaterial() {}

void TPZAxisymStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int64_t dimension = Dimension(); // problems dimension
    REAL radius = datavec[EVindex].x[0]; // radial distance to the axisymmetric z axis

    int64_t nShapeV = datavec[EVindex].fVecShapeIndex.NElements(); // number of velocity Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiV(dimension, nShapeV, 0.0);
    TPZFNMatrix<20, REAL>& divPhiV = datavec[EVindex].divphi;

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

    //dimension * (dimension + 1) / 2;
    TPZFNMatrix<150, REAL> StrainRate(3, nShapeV, 0.0); //Using voight notation
    TPZFNMatrix<150, REAL> StrainRateAxisymmetric(1, nShapeV, 0.0);
    
    if (datavec[EVindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {
            for (int64_t k = 0; k < 2; k++)
            {
                PhiV(k, j) = datavec[EVindex].fDeformedDirectionsFad(k, j).val();
            }

            StrainRate(0,j) = (datavec[EVindex].fDeformedDirectionsFad(0, j).fastAccessDx(0) - (PhiV(0,j) / radius)) / radius;
            StrainRate(1,j) = datavec[EVindex].fDeformedDirectionsFad(1, j).fastAccessDx(1) / radius;
            StrainRate(2,j) = 0.5 * (datavec[EVindex].fDeformedDirectionsFad(0, j).fastAccessDx(1) + datavec[EVindex].fDeformedDirectionsFad(1, j).fastAccessDx(0) - (PhiV(1,j) / radius)) / radius;
            StrainRateAxisymmetric(0,j) = PhiV(0,j) / (radius * radius);
        }
    }

    TPZFNMatrix<9, REAL> VoightCorrection(3, 3, 0.0); //This is to multiply by 2 the off diagonal part of strain rate tensor to account for its symmetry
    VoightCorrection(0,0) = 1.0;
    VoightCorrection(1,1) = 1.0;
    VoightCorrection(2,2) = 2.0;
    
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

    REAL factor = weight;
    ef.AddContribution(0, 0, PhiV, true, SourceTerm, false, factor);

    //Flux Matrix A contribution
    TPZFNMatrix<150, REAL> matrixC;
    VoightCorrection.Multiply(StrainRate, matrixC);
   
    factor = 2.0 * fviscosity * weight * radius;
    ek.AddContribution(0, 0, StrainRate, true, matrixC, false, factor);

    //Axisymmetric Flux Matrix A contribution
    ek.AddContribution(0, 0, StrainRateAxisymmetric, true, StrainRateAxisymmetric, false, factor);

    //Divergence Matrix B contribution
    factor = -1.0 * weight;
    ek.AddContribution(0, nShapeV, divPhiV, false, PhiP, true, factor);

    //Divergence Matrix BT contribution
    ek.AddContribution(nShapeV, 0, PhiP, false, divPhiV, true, factor);

    if(datavec.size()>2)
    {
        TPZFMatrix<REAL>& phivM = datavec[EVMindex].phi;
        TPZFMatrix<REAL>& phipM = datavec[EPMindex].phi;
        
        // Pressure and distributed flux
        for(int j=0; j<nShapeP; j++)
        {
            ek(nShapeV+nShapeP, nShapeV+j) += PhiP(j,0)*phivM(0,0)*weight;
            ek(nShapeV+j, nShapeV+nShapeP) += PhiP(j,0)*phivM(0,0)*weight;
        }
        
        // Injection and average-pressure
        ek(nShapeV+nShapeP+1, nShapeV+nShapeP) += phivM(0,0)*phipM(0,0)*weight;
        ek(nShapeV+nShapeP, nShapeV+nShapeP+1) += phivM(0,0)*phipM(0,0)*weight;
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

void TPZAxisymStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    TPZFNMatrix<150, REAL> PhiV = datavec[EVindex].phi;
    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    
    int64_t nShapeV = PhiV.Rows();
    int64_t nShapeP = PhiP.Rows();
            
    TPZFNMatrix<3,STATE> val1 = bc.Val1();
    TPZManVector<STATE, 3> val2 = bc.Val2();

    switch (bc.Type())
    {
    case 0: // Normal Velocity
    {
        REAL v_n = val2[0];
        REAL radius = datavec[EVindex].x[0];
        REAL factor = fBigNumber * weight;

        for (int64_t j = 0; j < nShapeV; j++)
        {
            ef(j) += v_n * PhiV(j,0) * factor;

            for (int64_t i = 0; i < nShapeV; i++)
            {
                ek(i, j) += PhiV(i,0) * PhiV(j,0) * factor / radius;
            }
        }
    }
    break;

    case 1: // Tangential Velocity
    {
        REAL v_t = val2[0];
        REAL radius = datavec[EPindex].x[0];
        REAL factor = weight * radius;

        for (int64_t i = 0; i < nShapeP; i++)
        {
            ef(i) += PhiP(i,0) * v_t * factor;
        }
    }
    break;

    case 2: // Normal Stress
    {
        REAL sigma_n = val2[0];

        for (int64_t i = 0; i < nShapeV; i++)
        {
            REAL phi = PhiV(i,0);
            ef(i) += sigma_n * PhiV(i,0) * weight;
        }
    }
    break;

    case 3: // Tangential Stress
    {
        REAL sigma_t = val2[0];
        REAL radius = datavec[EPindex].x[0];
        REAL factor = fBigNumber * weight * radius;

        for (int64_t j = 0; j < nShapeP; j++)
        {
            ef(j) += sigma_t * PhiP(j,0) * factor;

            for (int64_t i = 0; i < nShapeP; i++)
            {
                ek(i, j) += PhiP(i,0) * PhiP(j,0) * factor;
            }
        }
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

void TPZAxisymStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    TPZManVector<STATE, 3> v_h = datavec[EVindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];
    
    REAL radius = datavec[EVindex].x[0];
    if (abs(radius) < 1.0e-9) radius = 1.0e-6;
    
    Solout.Resize(NSolutionVariables(var));

    int64_t dimension = Dimension();
    
    switch(var) {
            
        case 0:{
            Solout[0] = p_h[0];
            break;
        }
            
        case 1:{
            Solout[0] = v_h[0] / radius; //Vx
            Solout[1] = v_h[1] / radius; //Vy
            Solout[2] = v_h[2] / radius; //Vz
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
            
                TPZFNMatrix<10,STATE> gradU = datavec[EVindex].dsol[0];
                gradU(2,2) = v_h[0] / radius;
                gradU *= 1.0/radius;

                gradU(0,0) -= (v_h[0] / (radius*radius));
                gradU(1,0) -= (v_h[1] / (radius*radius));
                
                TPZFNMatrix<9,STATE> StrainRate(3,3,0.), sigma(3,3,0.);
            
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        StrainRate(i,j)= 0.5 * (gradU(i,j) + gradU(j,i));
                        sigma(i,j) = 2.0 * fviscosity * StrainRate(i,j);
                    }
                    sigma(i,i) -= p_h[0];
                }

                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        Solout[i*3+j] = sigma(i,j);
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

int TPZAxisymStokesMaterial::IntegrationRuleOrder(const TPZVec<int>& elPMaxOrder) const
{
    const int maxOrder = [&elPMaxOrder=std::as_const(elPMaxOrder)](){
        int max = 0;
        for (auto ord : elPMaxOrder)
            if (ord > max) max = ord;
        return max;
    }();

    int ffporder = HasForcingFunction()? ForcingFunctionPOrder() : 0;

    const int intOrder = maxOrder < ffporder ? 10 * (maxOrder + ffporder) : (10 * maxOrder);

    return  intOrder;
}

int TPZAxisymStokesMaterial::IntegrationRuleOrderBC(const TPZVec<int>& elPMaxOrder) const
{
    const int maxOrder = [&elPMaxOrder=std::as_const(elPMaxOrder)](){
        int max = 0;
        for (auto ord : elPMaxOrder)
            if (ord > max) max = ord;
        return max;
    }();

    int ffporder = HasForcingFunction()? ForcingFunctionPOrder() : 0;

    const int intOrder = maxOrder < ffporder ? 10 * (maxOrder + ffporder) : (10 * maxOrder);

    return  intOrder;
}

void TPZAxisymStokesMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors){
    
    if(!HasExactSol()) DebugStop();

    errors.Resize(NEvalErrors());

    REAL radius = data[EVindex].x[0];
    
    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9,STATE> gradsol_exact(3,3);
    
    //Getting the exact solution for velocity, pressure and velocity gradient
    fExactSol(data[EVindex].x, sol_exact, gradsol_exact);
    
    //Getting the numeric solution for velocity, pressure and velocity gradient
    TPZManVector<STATE> v_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10,STATE> gradv_h = data[EVindex].dsol[0];
    
    this->Solution(data, VariableIndex("Velocity"), v_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);

    gradv_h(2, 2) = v_h[0] / radius;
    gradv_h *= 1.0 / radius;

    gradv_h(0, 0) -= (v_h[0] / (radius * radius));
    gradv_h(1, 0) -= (v_h[1] / (radius * radius));
    
    STATE diffv, diffp, diffdiv;

    diffp = p_h[0] - sol_exact[3];
    errors[0] = diffp * diffp;
    
    errors[1] = 0.0;
    for(int i = 0; i < 3; i++)
    {
        diffv = v_h[i] - sol_exact[i];
        errors[1] += diffv * diffv;
    }
    
    STATE div_exact = 0.0, div_h = 0.0;
    for(int i = 0; i < 3; i++)
    {
        div_exact += gradsol_exact(i,i);
        div_h += gradv_h(i, i);
    }
    
    diffdiv = div_h - div_exact;
    errors[2] = diffdiv * diffdiv;
}
