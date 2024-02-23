#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZMatInterfaceCombinedSpaces.h>
#include <TPZMatErrorCombinedSpaces.h>
#include <math.h>

#ifndef TPZSTOKESMATERIAL
#define TPZSTOKESMATERIAL

class TPZStokesMaterial : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>> {
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>>;
    
protected:
    /// material dimension
    int fdimension;
    
    /// Approximation Space for Velocity
    int fSpace = 1;
    
    /// fluid viscosity
    double fviscosity;
    
    enum SpaceIndex {EVindex, EPindex, EGindex, EPMindex};
    
//    /// Velocity index in datavec
//    int fVindex = 0;
//
//    /// Pressure index in datavec
//    int fPindex = 1;
//
//    int fPMindex = 2;
//
//    int fVMindex = 3;
    
    /// Big number for penalization method
    REAL fBigNumber = pow(10,std::numeric_limits<STATE>::max_digits10*2/3);
    
public:
    enum MVoigt {
        Exx, Exy, Eyy, Exz, Eyz, Ezz
    };
    
    /// Empty Constructor
    TPZStokesMaterial();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZStokesMaterial(int matID, int dimension, double viscosity);
    
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
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <class TVar>
    TVar TensorInnerProduct(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
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
    
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors) override;
    
    int NEvalErrors() const override {return 8;}
    
    virtual void StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6, REAL>& sigma, REAL pressure);
    
    virtual void StrainTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6, REAL>& epsilon);
    
    virtual void ViscosityTensor(TPZFNMatrix<36, REAL>& D);
    
    void ToVoigt(const TPZFNMatrix<9, STATE> &Sigma, TPZFNMatrix<6, STATE> &SVoigt) const;
    
    void FromVoigt(const TPZFNMatrix<6, STATE> &Svoigt, TPZFNMatrix<9, STATE> &S) const;
};

#endif
