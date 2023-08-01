#include "TPZ1dStokesMaterial.h"
#include <pzlog.h>

#ifdef PZ_LOG
static TPZLogger logger("pz.axisymmetricstokesmaterial");
#endif

TPZ1dStokesMaterial::TPZ1dStokesMaterial(): TBase() {}

TPZ1dStokesMaterial::TPZ1dStokesMaterial(int matID, REAL viscosity, REAL radius) : TBase(matID), fviscosity(viscosity), fradius(radius) {}

TPZ1dStokesMaterial::~TPZ1dStokesMaterial() {}

void TPZ1dStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{

#ifdef USING_LAPACK
    
    auto shapeType = datavec[EVindex].fShapeType;

    switch (datavec[EVindex].fShapeType)
    {
        case TPZShapeData::MShapeFunctionType::EVecandShape: //Hdiv
        case TPZShapeData::MShapeFunctionType::EVecShape:
        {
            int64_t nShapeV = datavec[EVindex].fVecShapeIndex.NElements(); // number of velocity Hdiv shape functions
            TPZFNMatrix<20, REAL>& divPhiV = datavec[EVindex].divphi;

            TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
            int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

            REAL factor = fradius * fradius * weight;
            
            //Divergence Matrix B contribution
            ek.AddContribution(0, nShapeV, divPhiV, false, PhiP, true, weight);

            //Divergence Matrix BT contribution
            ek.AddContribution(nShapeV, 0, PhiP, false, divPhiV, true, weight);
            break;
        }
        case TPZShapeData::MShapeFunctionType::EScalarShape: //H1
        {
            TPZFMatrix<REAL>& PhiV = datavec[EVindex].phi;
            int64_t nShapeV = PhiV.Rows(); // number of velocity H1 shape functions
            TPZFNMatrix<60, REAL>& GradPhiV = datavec[EVindex].dphix;

            TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
            int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

            REAL factor = fradius * fradius * weight;
            
            //Divergence Matrix B contribution
            ek.AddContribution(0, nShapeV, GradPhiV, true, PhiP, true, weight);

            //Divergence Matrix BT contribution
            ek.AddContribution(nShapeV, 0, PhiP, false, GradPhiV, false, weight);
            break;
        }
        default:
            break;
    }
    
#else

DebugStop();

#endif

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

void TPZ1dStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
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

        for (int64_t j = 0; j < nShapeV; j++)
        {
            ef(j) += v_n * fBigNumber;

            for (int64_t i = 0; i < nShapeV; i++)
            {
                ek(i, j) += fBigNumber;
            }
        }
    }
    break;
    case 2: // Normal Stress
    {
        REAL sigma_n = val2[0];
        REAL factor = fradius * fradius;

        ef(0) += sigma_n;
    }
    break;

    default:
    {
        std::cout << "ERROR: BOUNDARY NOT ALOUD FOR THIS ELEMENT" << std::endl;
        DebugStop();
    }
    break;
    }
}

int TPZ1dStokesMaterial::VariableIndex(const std::string& name) const {
    
    if(!strcmp("Pressure", name.c_str())) return 0;
    if(!strcmp("Velocity", name.c_str())) return 1;
    if(!strcmp("Force", name.c_str())) return 2;
    if(!strcmp("Tension", name.c_str())) return 3;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZ1dStokesMaterial::NSolutionVariables(int var) const{
    
    int aux;
    switch (var) {
        case 0: // pressure  [scalar]
        case 1: // velocity in the element direction
        case 2: // external force in the 
        case 3: // stress tensor
            aux = 1;
            break;
            
        default:
        std::cout << "\n\nVar index not implemented!!!\n\n";
        DebugStop();
    }
    return aux;
}

void TPZ1dStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {

    TPZManVector<STATE, 3> v_h = datavec[EVindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];
    
    Solout.Resize(NSolutionVariables(var));
    
    switch(var) {
            
        case 0:{
            Solout[0] = p_h[0];
            break;
        }
            
        case 1:{
            Solout[0] = v_h[1];
            break;
        }
            
        case 2:{
            TPZVec<STATE> f(3,0.);
            
            if(this->HasForcingFunction()){
                this->ForcingFunction()(datavec[EVindex].x, f);
            }
            
            Solout[0] = f[0];
            break;
        }
        case 3: {
            Solout[0] = p_h[0];
            break;
        }

        default:{
            std::cout << "\n\nVar index not implemented\n\n";
            DebugStop();
        }
    }
}