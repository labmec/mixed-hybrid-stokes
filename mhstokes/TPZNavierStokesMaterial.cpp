#include "TPZNavierStokesMaterial.h"

TPZNavierStokesMaterial::TPZNavierStokesMaterial() {}

TPZNavierStokesMaterial::TPZNavierStokesMaterial(int matID, int dimension, double viscosity) : TPZStokesMaterial(matID, dimension, viscosity) {}

TPZNavierStokesMaterial::~TPZNavierStokesMaterial() {}

void TPZNavierStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZStokesMaterial::Contribute(datavec, weight, ek, ef); //calls Contribute from Stokes Material

    //Add contribution of convective term (Oseen equation)
    REAL dim = Dimension();

    int64_t nShapeV = datavec[EVindex].fVecShapeIndex.NElements(); //number of velocity Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiV(dim, nShapeV, 0.0); //matrix to store velocity shape functions
    TPZFNMatrix<20, REAL>& divPhiV = datavec[EVindex].divphi;            

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi; //matrix to store pressure shape functions
    int64_t nShapeP = PhiP.Rows(); //number of pressure L2 shape functions

    const int nterms = dim*dim; //number of terms in the gradient matrix (using voight notation, that represents matrices in vector shape)
    
    if (datavec[EVindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {
            for (int64_t i = 0; i < dim; i++)
            {
                PhiV(i, j) = datavec[EVindex].fDeformedDirectionsFad(i, j).val(); // Directions on the deformed element. Used for H(div) and H(curl) approximation spaces.
            }
        }
    }

    //[[Exx, Exy, Eyx, Eyy, Exz, Eyz, Ezx, Ezy, Ezz], ...]^T;
    TPZFNMatrix<150, REAL> GradPhiV(nterms, nShapeV, 0.0); //Using voight notation
    if (datavec[EVindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {

            GradPhiV(Exx,j) = datavec[EVindex].fDeformedDirectionsFad(0, j).fastAccessDx(0); //fDeformedDirectionsFad(ui, j).fastAccessDx(xj) D(ui)/D(xj) --> uj = [u, v, w]  x = [x, y, z]
            GradPhiV(Exy,j) = datavec[EVindex].fDeformedDirectionsFad(0, j).fastAccessDx(1); 
            GradPhiV(Eyx,j) = datavec[EVindex].fDeformedDirectionsFad(1, j).fastAccessDx(0); 
            GradPhiV(Eyy,j) = datavec[EVindex].fDeformedDirectionsFad(1, j).fastAccessDx(1);
            if (dim > 2)
            {
                GradPhiV(Exz,j) = datavec[EVindex].fDeformedDirectionsFad(0, j).fastAccessDx(2);
                GradPhiV(Eyz,j) = datavec[EVindex].fDeformedDirectionsFad(1, j).fastAccessDx(2);
                GradPhiV(Ezx,j) = datavec[EVindex].fDeformedDirectionsFad(2, j).fastAccessDx(0);
                GradPhiV(Ezy,j) = datavec[EVindex].fDeformedDirectionsFad(2, j).fastAccessDx(1);
                GradPhiV(Ezz,j) = datavec[EVindex].fDeformedDirectionsFad(2, j).fastAccessDx(2); 
            }
        }
    }


    //TPZFNMatrix<150, REAL> gradPhiV(nterms, nShapeV, 0.0); //Using voight notation 
    //ToVoigt(GradPhiV, gradPhiV); //Converting to Voigt notation

    TPZManVector<REAL, 3> Beta(dim,0.); //Initializing Beta vector
    TPZFMatrix<REAL> gradBeta(dim,dim); 
    // if (this->HasExactSol()) // if an exact solution is given
    // {
    //     fExactSol(datavec[EVindex].x, Beta, gradBeta); //Assign Beta some exact solution
    // }
    // else
    // {
    //     if (fBeta) // if fBeta is set, call the function to compute Beta
    //     {
    //         fBeta(datavec[EVindex].x, Beta); //Assign Beta vector with the value of the beta function at the integration point
    //     }
    //     else
    //     {
    //         DebugStop();
    //     }
    // }
    if (fBeta) // if an exact solution is given
    {
        fBeta(datavec[EVindex].x, Beta); //Assign Beta vector with the value of the beta function at the integration point
    }
    else
    {
        if (this->HasExactSol()) // if fBeta is set, call the function to compute Beta
        {
            fExactSol(datavec[EVindex].x, Beta, gradBeta); //Assign Beta some exact solution
        }
        else
        {
            DebugStop();
        }
    }

    //[[Exx, Exy, Eyx], ...]^T;
    TPZFNMatrix<150, REAL> GradPhiVBeta(dim, nShapeV, 0.0);

    for (int64_t j = 0; j < nShapeV; j++)
    {

        GradPhiVBeta(Exx,j) = GradPhiV(Exx,j)*Beta[0] + GradPhiV(Exy,j)*Beta[1];
        GradPhiVBeta(Exy,j) = GradPhiV(Eyx,j)*Beta[0] + GradPhiV(Eyy,j)*Beta[1];

        if (dim > 2)
        {
            GradPhiVBeta(Exx,j) += GradPhiV(Exz,j)*Beta[2];
            GradPhiVBeta(Exy,j) += GradPhiV(Eyz,j)*Beta[2];
            GradPhiVBeta(Eyx,j) = GradPhiV(Ezx,j)*Beta[0] + GradPhiV(Ezy,j)*Beta[1] + GradPhiV(Ezz,j)*Beta[2];
        }
    }


    REAL factor = 1.0 * weight;
    ek.AddContribution(0, 0, PhiV, true, GradPhiVBeta, false, factor); // += weight * (PhiV^T.GradPhiVBeta))
}

    /** Receives a function that do some calculation on point x = [x, y, z]
     * and the result is stored in beta vector
     */
void TPZNavierStokesMaterial::setBetaFunc(TFunction betafunc){
    if (!betafunc) {
        DebugStop();
    }
    fBeta = betafunc;  
}

void TPZNavierStokesMaterial::ToVoigt(const TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt) const {

    Svoigt[Exx] = S(0, 0);
    Svoigt[Exy] = S(0, 1);
    Svoigt[Eyx] = S(1, 0);
    Svoigt[Eyy] = S(1, 1);
    if (Dimension() > 2) {
        Svoigt[Exz] = S(0, 2);
        Svoigt[Eyz] = S(1, 2);
        Svoigt[Ezx] = S(2, 0);
        Svoigt[Ezy] = S(2, 1);
        Svoigt[Ezz] = S(2, 2);
    }
}