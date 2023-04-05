#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include "TPZMatInterfaceCombinedSpaces.h"

#ifndef TPZSTOKESMATERIAL
#define TPZSTOKESMATERIAL

class TPZStokesMaterial : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceCombinedSpaces<STATE> > {
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceCombinedSpaces<STATE> >;
    
protected:
    /// material dimension
    int fdimension;
    
    /// Approximation Space for Velocity
    int fSpace = 1;
    
public:
    /// Empty Constructor
    TPZStokesMaterial();
    
    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZStokesMaterial(int matID, int dimension);
    
    /// Destructor
    ~TPZStokesMaterial();
    
    /// returns the solution associated with the var index based on the finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>>&datavec, int var, TPZVec<STATE>& Solout) override;
    
    /// returns the number of variables associated with the variable indexed by var. Var is obtained by calling VariableIndex
    int NSolutionVariables(int var) const override;
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    /** returns the integrable dimension of the material */
    int Dimension() const override {
        return fdimension;
    }
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const override {
        
        if(fSpace==1){
            return 1;
        }else{
            return 2;
        }
        
    } // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    
    virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                   const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                   int var, TPZVec<STATE> &Solout) override {
        DebugStop();
    }
    
    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation around
     * one interface element
     */
    virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                   const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                   int var, TPZVec<STATE> &Solout,
                                   TPZCompEl *left,TPZCompEl *right) override
    {
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * */
    virtual void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                          const std::map<int, TPZMaterialDataT<STATE>> &datavec,
                          REAL weight,
                          TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                                   TPZBndCondT<STATE> &bc) override;
    
    /**
     * @brief Computes a contribution to the stiffness matrix and load vector at one integration point
     * to multiphysics simulation
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ek is the stiffness matrix
     * @param [out] ef is the load vector
     * */
    virtual void ContributeInterface(const TPZMaterialDataT<STATE> &data,
                                     const std::map<int, TPZMaterialDataT<STATE>> &dataleft,
                                     const std::map<int, TPZMaterialDataT<STATE>> &dataright, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef) override;
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution
     * of interface elements
     */
    virtual void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                  std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                  std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;
};

#endif
