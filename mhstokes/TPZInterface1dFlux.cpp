#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzlog.h>

#include "TPZInterface1dFlux.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesInterface");
#endif


TPZInterface1dFlux::TPZInterface1dFlux(int matID, REAL viscosity, REAL radius) : TBase(matID), fviscosity(viscosity), fradius(radius) {}

TPZInterface1dFlux::~TPZInterface1dFlux() {}

void TPZInterface1dFlux::ContributeInterface(const TPZMaterialDataT<STATE>& data, const std::map<int, TPZMaterialDataT<STATE>>& dataleft, const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef)
{
    if(dataleft.find(fVindex) == dataleft.end()) DebugStop();
    if(dataright.find(fPindex) == dataright.end()) DebugStop();
    
    const TPZMaterialDataT<STATE>& vDataLeft = dataleft.find(fVindex)->second;
    const TPZMaterialDataT<STATE>& pDataRight = dataright.find(fPindex)->second;

    switch (vDataLeft.fShapeType)
    {
        case TPZShapeData::MShapeFunctionType::EVecandShape: //Hdiv
        case TPZShapeData::MShapeFunctionType::EVecShape:
        {
            int64_t nShapeV = vDataLeft.fVecShapeIndex.NElements(); // number of velocity Hdiv shape functions
            TPZFNMatrix<3, REAL> phiV(nShapeV, 1, 0.0);
            
            for (int64_t i = 0; i < nShapeV; i++)
            {
                phiV(i,0) = vDataLeft.fDeformedDirections(0,i); //only the x component
            }

            TPZFNMatrix<3, REAL> phiP = pDataRight.phi;
            int64_t nShapePressure = pDataRight.phi.Rows(); // number of pressure H1 shape functions

            REAL factor = fMultiplier * weight;

            ek.AddContribution(0, nShapeV, phiV, false, phiP, true, factor);

            ek.AddContribution(nShapeV, 0, phiP, false, phiV, true, factor);
            break;
        }
        case TPZShapeData::MShapeFunctionType::EScalarShape: //H1
        {
            int64_t nShapeV = vDataLeft.phi.Rows();
            int64_t nShapeLambda = pDataRight.phi.Rows();
            
            TPZFNMatrix<3, REAL> phiV = vDataLeft.phi;
            TPZFNMatrix<3, REAL> phiP = pDataRight.phi;

            REAL factor = fMultiplier * weight;

            ek.AddContribution(0, nShapeV, phiV, false, phiP, true, factor);

            ek.AddContribution(nShapeV, 0, phiP, false, phiV, true, factor);
            break;
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

void TPZInterface1dFlux::ContributeBCInterface(const TPZMaterialDataT<STATE> &data, const std::map<int, TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{    
    DebugStop();
}

void TPZInterface1dFlux::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZInterface1dFlux::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  REAL weight, TPZFMatrix<STATE> &ek,
                  TPZFMatrix<STATE> &ef,
                  TPZBndCondT<STATE> &bc)
{
    DebugStop();
}

void TPZInterface1dFlux::FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, std::map<int, TPZMaterialDataT<STATE>> &datavec_right)
{
    datavec_left[0].fNeedsNormal = true;
    datavec_left[0].fNeedsSol = true;
    datavec_right[1].fNeedsNormal = true;
    datavec_right[1].fNeedsSol = true;
    datavec_left[0].fNeedsDeformedDirectionsFad = false;
}

void TPZInterface1dFlux::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                                     const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                                     const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                                     int var, TPZVec<STATE> &Solout)
{
    DebugStop();
}

void TPZInterface1dFlux::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                                     const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                                     const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                                     int var, TPZVec<STATE> &Solout,
                                                     TPZCompEl *left, TPZCompEl *right)
{
    DebugStop();
}

void TPZInterface1dFlux::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &sol)
{

    DebugStop();
}
