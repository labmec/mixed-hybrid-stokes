#include "TPZStokesMaterial.h"

#ifndef TPZAXISYMSTOKESMATERIAL
#define TPZAXISYMSTOKESMATERIAL

class TPZAxisymStokesMaterial : public TPZStokesMaterial {
    
public:
    /// Empty Constructor
    TPZAxisymStokesMaterial();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZAxisymStokesMaterial(int matID, int dimension, double viscosity);
    
    /// Destructor
    ~TPZAxisymStokesMaterial();
    
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

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) override;

    [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;

    [[nodiscard]] int IntegrationRuleOrderBC(const TPZVec<int> &elPMaxOrder) const override;

    virtual int NEvalErrors() const override {return 3;}

    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors) override;

};

#endif