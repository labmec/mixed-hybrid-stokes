#include <pzfmatrix.h>
#include <TPZBndCondT.h>

#include "TPZStokesMaterial.h"

TPZStokesMaterial::TPZStokesMaterial(): TBase() {
    
}

TPZStokesMaterial::TPZStokesMaterial(int matID, int dimension) : TBase(matID), fdimension(dimension){
    
}

TPZStokesMaterial::~TPZStokesMaterial(){
    
}

void TPZStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef){
    DebugStop();
    // YOU NEED TO IMPLEMENT THIS
}

void TPZStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    DebugStop();
    // AND THIS
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
    
    int Vindex = 0;
    int Pindex = 1;
    
    TPZManVector<STATE, 3> v_h = datavec[Vindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[Pindex].sol[0];
    
    TPZFNMatrix<9, STATE> gradU(3,1);
    
    TPZFMatrix<STATE>& dsol = datavec[Vindex].dsol[0];
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
                this->ForcingFunction()(datavec[Vindex].x, f);
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
    DebugStop();
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
