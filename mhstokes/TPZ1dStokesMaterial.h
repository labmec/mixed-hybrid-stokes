#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <math.h>

#ifndef TPZ1DSTOKESMATERIAL
#define TPZ1DSTOKESMATERIAL

class TPZ1dStokesMaterial : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>>
{
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>>;

    enum SpaceIndex{EVindex, EPindex};

protected:
    REAL fviscosity;
    REAL fradius;
    
public:
    /// Empty Constructor
    TPZ1dStokesMaterial();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZ1dStokesMaterial(int matID, REAL viscosity, REAL radius);
    
    /// Destructor
    ~TPZ1dStokesMaterial();
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    virtual int VariableIndex(const std::string& name) const override;

    virtual int NSolutionVariables(int var) const override;

    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) override;

    [[nodiscard]] int Dimension() const override {return 1;}

    [[nodiscard]] int NStateVariables() const override {return 1;}
};

#endif