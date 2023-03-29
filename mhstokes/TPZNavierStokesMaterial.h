///*
// *  TPZNavierStokesMaterial.cpp
// *  PZ
// *
// *  Created by Pablo Carvalho on 10/05/2016.
// *  Copyright 2016 __MyCompanyName__. All rights reserved.
// *
// */
//
//#include "TPZMatWithMem.h"
//#include "pzfmatrix.h"
//#include "TPZBndCondT.h"
//#include "pzlog.h"
//#include "tpzautopointer.h"
//#include "TPZMaterial.h"
//#include "pztrnsform.h"
//#include "TPZAnalyticSolution.h"
//#include "TPZSimulationData.h"
//#include "TPZNSMemory.h"
//
//#include "TPZMatBase.h"
//#include "TPZMatCombinedSpaces.h"
//#include "TPZMatErrorCombinedSpaces.h"
//#include "TPZMatInterfaceCombinedSpaces.h"
//#include "TPZMaterialDataT.h"
//#include "TPZMatWithMem.h"
//
//
//#ifndef TPZNavierStokesMATERIAL
//#define TPZNavierStokesMATERIAL
//
////enum NSProblemType {ENavierStokes,EOseen,EStokes,EBrinkman};
//
//class TPZNavierStokesMaterial : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TPZNSMemory>, TPZMatErrorCombinedSpaces<STATE>,TPZMatInterfaceCombinedSpaces<STATE> >{
//
//    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatWithMem<TPZNSMemory>, TPZMatErrorCombinedSpaces<STATE>, TPZMatInterfaceCombinedSpaces<STATE> >;
//
//protected:
//
//    /// dimension of the material
//    int fDimension;
//
//    /// Aproximation Space for velocity
//    int fSpace;
//
//    /// viscosidade
//    STATE fViscosity;
//
//    /// Brinkman coef
//    STATE fcBrinkman;
//
//    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
//    STATE fk;
//
//    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
//    STATE fTheta;
//
//    STATE fSigma;
//
//    TStokesAnalytic::MProblemType f_problemtype;
//
//    /** @brief Simulation time step */
//    REAL fDeltaT;
//
//    /** @brief State: one ou one+1 */
//    enum EState { ELastState = 0, ECurrentState = 1 };
//
//    EState fState;
//
//    /** Data for simulation */
//    TPZSimulationData *f_sim_data;
//
//public:
//
//    bool NeedsNormalVecFad = true;
//    /**
//     * Empty Constructor
//     */
//    TPZNavierStokesMaterial();
//
//    /** Creates a material object and inserts it in the vector of
//     *  material pointers of the mesh.
//     */
//    TPZNavierStokesMaterial(int matid, int dimension);
//
//
//    /** Creates a material object based on the referred object and
//     *  inserts it in the vector of material pointers of the mesh.
//     */
//    TPZNavierStokesMaterial(const TPZNavierStokesMaterial &mat);
//
//    /**
//     * Destructor
//     */
//    ~TPZNavierStokesMaterial();
//
//    TPZNavierStokesMaterial &operator=(const TPZNavierStokesMaterial &copy)
//    {
//        DebugStop();
//        return *this;
//    }
//
//    /** Fill material data parameter with necessary requirements for the
//     * Contribute method. Here, in base class, all requirements are considered
//     * as necessary. Each derived class may optimize performance by selecting
//     * only the necessary data.
//     * @since April 10, 2007
//     */
//    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
//
//    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
//    virtual void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
//
////    void SetTranform(TPZTransform<STATE> Transf, TPZTransform<STATE> InvTransf)
////    {
////        f_T = Transf;
////        f_InvT = InvTransf;
////    }
//
//    void SetSimulationData(TPZSimulationData *simdata);
//
//    void SetLastState(){ fState = ELastState; }
//
//    void SetCurrentState(){ fState = ECurrentState; }
//
//    void SetTimeStep(REAL timeStep) {
//        fDeltaT = timeStep;
//    }
//
//    void SetPermeability(REAL perm) {
//        fk = perm;
//    }
//
//    void SetProblemType(TStokesAnalytic::MProblemType type){
//        f_problemtype = type;
//    }
//
//    TStokesAnalytic::MProblemType GetProblemType(){
//        return f_problemtype;
//    }
//
//    /** returns the name of the material */
//    std::string Name() {
//        return "TPZNavierStokesMaterial";
//    }
//
//    STATE GetViscosity() const{
//        return fViscosity;
//    }
//
//    void SetViscosity(STATE visco){
//        fViscosity = visco;
//    }
//
//    STATE GetBrinkman() const{
//        return fcBrinkman;
//    }
//
//    void SetBrinkman(STATE cBrinkman){
//        fcBrinkman = cBrinkman;
//    }
//
//    /** returns the integrable dimension of the material */
//    int Dimension() const override {
//        return fDimension;
//    }
//
//    /** returns the number of state variables associated with the material */
//    virtual int NStateVariables() const override {
//
//        if(fSpace==1){
//            return 1;
//        }else{
//            return 2;
//        }
//
//    } // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
//
//    /** print out the data associated with the material */
//    void Print(std::ostream &out = std::cout);
//
//    /** returns the variable index associated with the name */
//    int VariableIndex(const std::string &name) const override;
//
//    /** returns the number of variables associated with the variable
//     indexed by var.  var is obtained by calling VariableIndex */
//    int NSolutionVariables(int var) const override;
//
//    /** Computes the divergence over the parametric space */
//    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialDataT<STATE>> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
//
//
//    /** returns the solution associated with the var index based on
//     * the finite element approximation */
//    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;
//
//    /** index of velocity */
//    int VIndex(){ return 0; }
//
//    /** index of pressure */
//    int PIndex(){ return 1; }
//
//    /** inner product of two tensors. See Gurtin (2003), p. 5. */
//    template <class TVar>
//    TVar Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
//
//    /** inner product of two vectors. See Gurtin (2003), p. 5. */
//    STATE InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);
//
//    /** trace of the tensor GradU = Div(U)*/
//    STATE Tr(TPZFMatrix<STATE> &GradU );
//
//    /** transpose of the tensor GradU = Div(U)*/
//    STATE Transpose(TPZFMatrix<STATE> &GradU );
//
//    /** Fill the vector of gradient for each phi */
//    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
//
//    /// transform a H1 data structure to a vector data structure
//    void FillVecShapeIndex(TPZMaterialData &data);
//
//
//    // Contribute Methods being used - Multiphysics
//    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
//                            REAL weight,TPZFMatrix<STATE> &ek,
//                            TPZFMatrix<STATE> &ef) override;
////    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
//
//
//    // Contribute Methods being used
//
//    /**
//     * It computes a contribution to the stiffness matrix and load vector at one integration point.
//     * @param data[in] stores all input data
//     * @param weight[in] is the weight of the integration rule
//     * @param ek[out] is the stiffness matrix
//     * @param ef[out] is the load vector
//     * @since April 16, 2007
//     */
//
//    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
//
//    /**
//     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
//     * @param data[in] stores all input data
//     * @param weight[in] is the weight of the integration rule
//     * @param ek[out] is the stiffness matrix
//     * @param ef[out] is the load vector
//     * @param bc[in] is the boundary condition material
//     * @since April 16, 2007
//     */
//    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
//
//    /**
//     * Save the element data to a stream
//     */
//    virtual void Write(TPZStream &buf, int withclassid) const override;
//
//    /**
//     * Read the element data from a stream
//     */
//    void Read(TPZStream &buf, void *context) override;
//
//
//    virtual int NEvalErrors() {return 6;}
//
//    /**
//     * It computes errors contribution in differents spaces.
//     * @param data[in] stores all input data
//     * @param weight[in] is the weight of the integration rule
//     * @param ef[out] is the load vector
//     * @param bc[in] is the boundary condition material
//     */
//    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
//
//};
//
//#endif
