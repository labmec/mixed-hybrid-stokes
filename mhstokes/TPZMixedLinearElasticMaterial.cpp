#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>

#include "TPZMixedLinearElasticMaterial.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.elasticmaterial");
#endif

TPZMixedLinearElasticMaterial::TPZMixedLinearElasticMaterial(): TBase() {}

TPZMixedLinearElasticMaterial::TPZMixedLinearElasticMaterial(int matID, int dimension, REAL young_modulus, REAL poisson, AnalysisType analysisType) : TBase(matID), fdimension(dimension), fyoung(young_modulus), fpoisson(poisson), fAnalysisType(analysisType) {}

TPZMixedLinearElasticMaterial::~TPZMixedLinearElasticMaterial() {}


void TPZMixedLinearElasticMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef)
{
    int64_t nShapeU = datavec[EUindex].fVecShapeIndex.NElements(); // number of displacements Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiU(fdimension, nShapeU, 0.0);
    TPZFNMatrix<20, REAL>& divPhiU = datavec[EUindex].divphi;

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure L2 shape functions

    const int nterms = fdimension * (fdimension + 1) / 2;
    TPZFNMatrix<150, REAL> Strain(nterms, nShapeU, 0.0); //Using voight notation
    
    if (datavec[EUindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeU; j++)
        {
            int cont = fdimension-1;
            for (int i = 0; i < fdimension; i++)
            {
                PhiU(i, j) = datavec[EUindex].fDeformedDirectionsFad(i, j).val();
                Strain(i,j) = datavec[EUindex].fDeformedDirectionsFad(i, j).fastAccessDx(i); //diagonal part of infinitesimal strain tensor
                for (int64_t k = i+1; k < fdimension && k != i; k++)
                {
                    Strain(++cont,j) = 0.5*(datavec[EUindex].fDeformedDirectionsFad(i, j).fastAccessDx(k) + datavec[EUindex].fDeformedDirectionsFad(k, j).fastAccessDx(i)); //off diagonal parte of infinitesimal strain tensor
                }
            }
        }
    }

    const REAL lame = fyoung*fpoisson / ((1.0+fpoisson)*(1.0-2.0*fpoisson));
    const REAL mu = 0.5 * fyoung / (1.0+fpoisson);
    const REAL bulk = lame + 2.0 * mu / 3.0;

    TPZFNMatrix<36, REAL> elasticityD(nterms, nterms, 0.0); // The terms related to the off diagonal part of strain tensor is multiplied by 2 to account for its symmetry

    switch (fAnalysisType)
    {
        case AnalysisType::EGeneral:
        case AnalysisType::EPlaneStrain:
        {
            for (int64_t i = 0; i < fdimension; i++)
            {
                elasticityD(i, i) = 2.0/3.0;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    elasticityD(i, j) = -1.0 / 3.0; //volumetric part
                    elasticityD(j, i) = -1.0 / 3.0;
                }
            }
            for (int64_t k = fdimension; k < nterms; k++)
            {
                elasticityD(k, k) = 2.0; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
        case AnalysisType::EPlaneStress:
        {
            for (int64_t i = 0; i < fdimension; i++)
            {
                elasticityD(i, i) = 2.0/3.0;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    elasticityD(i, j) = -1.0 / 3.0; //volumetric part
                    elasticityD(j, i) = -1.0 / 3.0;
                }
            }
            for (int64_t k = fdimension; k < nterms; k++)
            {
                elasticityD(k, k) = 2.0; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
    }
    
    TPZFNMatrix<3,REAL> SourceTerm(fdimension, 1.0, 0.0);
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EUindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
        {
            SourceTerm(i,0) = sourceAux[i];
        }
    }

    //Body Forces contribution
    ef.AddContribution(0, 0, PhiU, true, SourceTerm, false, weight);

    //Stiffness Matrix A contribution
    TPZFNMatrix<150, REAL> aux;
    elasticityD.Multiply(Strain, aux);
   
    REAL factor = 2.0 * mu * weight;
    ek.AddContribution(0, 0, Strain, true, aux, false, factor);

    //Divergence Matrix B contribution
    factor = -1.0 * weight;
    ek.AddContribution(0, nShapeU, divPhiU, false, PhiP, true, factor);

    //Divergence Matrix BT contribution
    ek.AddContribution(nShapeU, 0, PhiP, false, divPhiU, true, factor);

    //Bulk Matrix C contribution
    factor = -1.0 / bulk * weight;
    ek.AddContribution(nShapeU, nShapeU, PhiP, false, PhiP, true, factor);
    
    if(datavec.size() > 2) //Static condensation
    {
        TPZFMatrix<REAL>& PhiUM = datavec[EVMindex].phi;
        TPZFMatrix<REAL>& phipM = datavec[EPMindex].phi;
        
        // Pressure and distributed displacement
        for(int j = 0; j < nShapeP; j++)
        {
            ek(nShapeU + nShapeP, nShapeU + j      ) += PhiP(j,0) * PhiUM(0,0) * weight;
            ek(nShapeU + j      , nShapeU + nShapeP) += PhiP(j,0) * PhiUM(0,0) * weight;
        }
        
        // Injection and average-pressure
        ek(nShapeU + nShapeP + 1, nShapeU + nShapeP    ) += PhiUM(0,0) * phipM(0,0) * weight;
        ek(nShapeU + nShapeP    , nShapeU + nShapeP + 1) += PhiUM(0,0) * phipM(0,0) * weight;
    }

#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        elasticityD.Print("tensor D", sout, EMathematicaInput);
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZMixedLinearElasticMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    TPZFNMatrix<150, REAL> PhiU = datavec[EUindex].phi;
    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    
    int64_t nShapeU = PhiU.Rows();
    int64_t nShapeP = PhiP.Rows();
            
    TPZFNMatrix<3,STATE> val1 = bc.Val1();
    TPZManVector<STATE, 3> val2 = bc.Val2();

    switch (bc.Type())
    {
    case 0: // Normal Displacement
    {
        REAL u_n = val2[0];
        REAL factor = fBigNumber * weight;

        for (int64_t j = 0; j < nShapeU; j++)
        {
            ef(j) += u_n * PhiU(j,0) * factor;

            for (int64_t i = 0; i < nShapeU; i++)
            {
                ek(i, j) += PhiU(i,0) * PhiU(j,0) * factor;
            }
        }
    }
    break;

    case 1: // Tangential displacement
    {
        REAL u_t = val2[0];

        for (int64_t i = 0; i < nShapeP; i++)
        {
            ef(i) += PhiP(i,0) * u_t * weight;
        }
    }
    break;

    case 2: // Normal Stress
    {
        REAL sigma_n = val2[0];

        for (int64_t i = 0; i < nShapeU; i++)
        {
            REAL phi = PhiU(i,0);
            ef(i) += sigma_n * PhiU(i,0) * weight;
        }
    }
    break;

    case 3: // Tangential Stress
    {
        REAL sigma_t = val2[0];
        REAL factor = fBigNumber * weight;

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

int TPZMixedLinearElasticMaterial::VariableIndex(const std::string& name) const {
    
    if(!strcmp("Pressure", name.c_str())) return 0;
    if(!strcmp("Displacement", name.c_str())) return 1;
    if(!strcmp("Force", name.c_str())) return 2;
    if(!strcmp("Stress", name.c_str())) return 3;
    if(!strcmp("Strain", name.c_str())) return 4;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZMixedLinearElasticMaterial::NSolutionVariables(int var) const{
    
    int aux;
    switch (var) {
        case 0: // pressure  [scalar]
            aux = 1;
            break;
        case 1: // displacement [vector]
        case 2: // external force [vector]
            aux = 3;
            break;
        case 3: // stress tensor
        case 4: // strain tensor
            aux = 9;
            break;
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
            break;
    }
    return aux;
}

void TPZMixedLinearElasticMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    TPZManVector<STATE, 3> u_h = datavec[EUindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];

    const REAL lame = fyoung*fpoisson / ((1.0+fpoisson)*(1.0-2.0*fpoisson));
    const REAL mu = 0.5 * fyoung / (1.0+fpoisson);
    
    Solout.Resize(NSolutionVariables(var));
    
    switch(var) {
            
        case 0:
        {
            Solout[0] = p_h[0];
            break;
        }
            
        case 1:
        {
            Solout[0] = u_h[0]; //Vx
            Solout[1] = u_h[1]; //Vy
            Solout[2] = u_h[2]; //Vz
            break;
        }
            
        case 2:
        {
            TPZVec<STATE> f(3,0.);
            
            if(this->HasForcingFunction()){
                this->ForcingFunction()(datavec[EUindex].x, f);
            }
            
            Solout[0] = f[0];
            Solout[1] = f[1];
            Solout[2] = f[2];
            break;
        }
        case 3:
        {

            TPZFNMatrix<10, STATE> gradU = datavec[EUindex].dsol[0];
            TPZFNMatrix<9, STATE> strain(3, 3, 0.0), sigma(3, 3, 0.0);

            REAL trace = 0;
            for (int i = 0; i < fdimension; i++)
            {
                for (int j = 0; j < fdimension; j++)
                {
                    strain(i, j) = 0.5 * (gradU(i, j) + gradU(j, i));
                }
                trace += strain(i,i);
            }

            for (int i = 0; i < fdimension; i++)
            {
                for (int j = 0; j < fdimension; j++)
                {
                    sigma(i,j) = 2.0 * mu * strain(i,j);
                }
                sigma(i,i) -= 2.0 / 3.0 * mu * trace + p_h[0];
            }
            
            if (fAnalysisType == AnalysisType::EPlaneStrain)
                sigma(2,2) = fpoisson * (sigma(0,0) + sigma(1,1));

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Solout[i * 3 + j] = sigma(i, j);
                }
            }
            break;
        }
        case 4:
        {
            TPZFNMatrix<10, STATE> gradU = datavec[EUindex].dsol[0];
            TPZFNMatrix<9, STATE> strain(3, 3, 0.0), sigma(3, 3, 0.0);

            for (int i = 0; i < fdimension; i++)
            {
                for (int j = 0; j < fdimension; j++)
                {
                    strain(i, j) = 0.5 * (gradU(i, j) + gradU(j, i));
                }
            }

            if (fAnalysisType == AnalysisType::EPlaneStress)
                strain(2, 2) = -lame * (strain(0,0) + strain(1,1)) / (lame + 2.0 * mu);
        }
 
        default:{
            std::cout << "\n\nVar index not implemented\n\n";
            DebugStop();
        }
    }
}

void TPZMixedLinearElasticMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
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

void TPZMixedLinearElasticMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors){
    
    if(!HasExactSol()) DebugStop();
    
    TPZManVector<STATE, 4> sol_exact(3);
    TPZFNMatrix<9> dsol_exact(3,3);
    
    fExactSol(data[EUindex].x, sol_exact, dsol_exact);
    
    errors.Resize(NEvalErrors());
    errors.Fill(0.);
    
    TPZManVector<STATE> Displacement(3,0.);
    TPZManVector<STATE> Pressure(3,0.);
    
    this->Solution(data, VariableIndex("Displacement"), Displacement);
    this->Solution(data, VariableIndex("Pressure"), Pressure);
    
    TPZFMatrix<STATE>& dsolv = data[EUindex].dsol[0];
    
    dsolv.Resize(3, 3);
    
    STATE diffv, diffp, diffdiv;
    
    errors[0] = 0.;
    for(int i=0; i<3; i++){
        diffv = Displacement[i] - sol_exact[i];
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
