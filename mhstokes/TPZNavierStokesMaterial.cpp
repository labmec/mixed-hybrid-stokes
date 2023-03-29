///*
// *  TPZNavierStokesMaterial.cpp
// *  PZ
// *
// *  Created by Pablo Carvalho on 10/05/2016.
// *  Copyright 2016 __MyCompanyName__. All rights reserved.
// *
// */
//
//#include "TPZNavierStokesMaterial.h"
//#include "TPZBndCondT.h"
//#include "pzaxestools.h"
//#include "TPZMatWithMem.h"
//#include "pzfmatrix.h"
//#include "pzlog.h"
//
//using namespace std;
//
//
//TPZNavierStokesMaterial::TPZNavierStokesMaterial() : TBase() {
//    //fDim = 1;
//    TPZFNMatrix<3,STATE> Vl(1,1,0.);
//    TPZNSMemory defaultmemory;
//    this->SetDefaultMem(defaultmemory);
//    fk=1;
//    fViscosity=1.;
//    fcBrinkman=0.;
//    f_problemtype = TStokesAnalytic::ENavierStokes;
//    fState = ECurrentState;
//    fDeltaT = 0.;
//
//}
//
//////////////////////////////////////////////////////////////////////
//
//TPZNavierStokesMaterial::TPZNavierStokesMaterial(int matid, int dimension) : TBase(matid),fDimension(dimension),fSpace(1),fViscosity(1.),fcBrinkman(0.),fTheta(0),fSigma(0)
//{
//    // symmetric version
//    //fTheta = -1;
//
//    TPZFNMatrix<3,STATE> Vl(1,1,0.);
//    TPZNSMemory defaultmemory;
//    this->SetDefaultMem(defaultmemory);
//    fk=1.;
//    f_problemtype = TStokesAnalytic::ENavierStokes;
//    fState = ECurrentState;
//    fDeltaT = 0.;
//    
//}
//
//////////////////////////////////////////////////////////////////////
//
//TPZNavierStokesMaterial::TPZNavierStokesMaterial(const TPZNavierStokesMaterial &mat) : TBase(mat),fDimension(mat.fDimension),fSpace(mat.fSpace), fViscosity(mat.fViscosity),fcBrinkman(mat.fcBrinkman), fTheta(mat.fTheta), fSigma(mat.fSigma)
//{
//    fk= mat.fk;
//    f_problemtype = mat.f_problemtype;
//    fState = mat.fState;
//    fDeltaT = mat.fDeltaT;
//}
//
//////////////////////////////////////////////////////////////////////
//
//TPZNavierStokesMaterial::~TPZNavierStokesMaterial(){
//    
//    
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::SetSimulationData(TPZSimulationData *simdata){
//    f_sim_data = simdata;
//
//    //Set parameters and constants:
//    f_problemtype = f_sim_data->GetProblemType();
//    fViscosity = f_sim_data->GetViscosity();
//    fcBrinkman = f_sim_data->GetBrinkmanCoef();
//}
//
//////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
//{
//    int ndata = datavec.size();
//    for (int idata=0; idata < ndata ; idata++) {
//        datavec[idata].SetAllRequirements(false);
//        datavec[idata].fNeedsSol = true;
//        datavec[idata].fNeedsHSize = true;
//        datavec[idata].fNeedsNormal = true;
//    }
//    datavec[0].fNeedsDeformedDirectionsFad = NeedsNormalVecFad;
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE>> &datavec) const
//{
//    int ndata = datavec.size();
//    for (int idata=0; idata < ndata ; idata++) {
//        datavec[idata].SetAllRequirements(false);
//        datavec[idata].fNeedsSol = true;
//        datavec[idata].fNeedsNormal = true;
//        datavec[idata].fNeedsNeighborSol = true;
//    }
//    datavec[0].fNeedsDeformedDirectionsFad = NeedsNormalVecFad;
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::Print(std::ostream &out) {
//    out << "\t Base class print:\n";
//    out << " name of material : " << this->Name() << "\n";
//    TPZMaterial::Print(out);
//}
//
//////////////////////////////////////////////////////////////////////
//
//int TPZNavierStokesMaterial::VariableIndex(const std::string &name) const {
//    
//    if (!strcmp("P", name.c_str()))  return 0;
//    if (!strcmp("Pressure", name.c_str()))  return 0;
//    if (!strcmp("V", name.c_str()))  return 1;
//    if (!strcmp("State", name.c_str()))  return 0;
//    if (!strcmp("f", name.c_str()))         return 2;
//    if (!strcmp("V_exact", name.c_str()))   return 3;
//    if (!strcmp("P_exact", name.c_str()))   return 4;
//    if (!strcmp("Div", name.c_str()))   return 5;
//    if (!strcmp("SymTensorNorm", name.c_str()))   return 6;
//    if (!strcmp("P_CDG", name.c_str()))  return 7;
//    if (!strcmp("P_exact_CDG", name.c_str()))  return 8;
//    if (!strcmp("Vorticity2D", name.c_str()))  return 9;
//    if (!strcmp("Tension",name.c_str())) return 10;
//    //    if (!strcmp("V_exactBC", name.c_str()))   return 5;
//    
//    std::cout  << " Var index not implemented " << std::endl;
//    DebugStop();
//    return 0;
//}
//
//////////////////////////////////////////////////////////////////////
//
//int TPZNavierStokesMaterial::NSolutionVariables(int var) const {
//    
//    switch(var) {
//            
//        case 0:
//            return 1; // Pressure, Scalar
//        case 1:
//            return 3; // Velocity, Vector
//        case 2:
//            return 3; // f, Vector
//        case 3:
//            return 3; // V_exact, Vector
//        case 4:
//            return 1; // P_exact, Scalar
//        case 5:
//            return 1; // Divergente
//        case 6:
//            return 1; // Symetric tensor norm
//        case 7:
//            return 1; // Pressure for CDG Navier-Stokes formulation
//        case 8:
//            return 1; // Exact pressure for CDG Navier-Stokes formulation
//        case 9:
//            return 1; // 2D - Vorticity (z-direction)
//        case 10:
//            return 9;
//
//            //        case 5:
//            //            return this->Dimension(); // V_exactBC, Vector
//        default:
//        {
//            std::cout  << " Var index not implemented " << std::endl;
//            DebugStop();
//        }
//    }
//    return 0;
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {
//    
//    
//    int vindex = this->VIndex();
//    int pindex = this->PIndex();
//    
//    TPZManVector<STATE,3> v_h = datavec[vindex].sol[0];
//    TPZManVector<STATE,3> p_h = datavec[pindex].sol[0];
//    
//    TPZFNMatrix<9,STATE> gradu(3,1);
//    
//    // TPZManVector<STATE> v_h = datavec[vindex].sol[0];
//    // TPZManVector<STATE> p_h = datavec[pindex].sol[0];
//    
//    TPZFMatrix<STATE> &dsol = datavec[vindex].dsol[0];
//   // dsol.Resize(3,3);
//    TPZFNMatrix<9,STATE> dsolxy(3,3),dsolxyp(3,1);
//    dsolxy = dsol;
//    if (fSpace!=1) {
//        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, datavec[vindex].axes);
//    }
//
//    Solout.Resize(this->NSolutionVariables(var));
//    
//    switch(var) {
//            
//        case 0: //Pressure
//        {
//            Solout[0] = p_h[0];
//        }
//            break;
//
//        case 1: //Velocity
//        {
//            Solout[0] = v_h[0]; // Vx
//            Solout[1] = v_h[1]; // Vy
//            Solout[2] = v_h[2]; // Vz
//        }
//            break;
//        case 2: //f
//        {
//            TPZVec<STATE> f(3,0.0);
//            if(f_problemtype==TStokesAnalytic::EBrinkman){
//                f.resize(4);
//            }
//            if(this->HasForcingFunction()){
//                TPZVec<STATE> x(3,0.);
//                x=datavec[vindex].x;
//                this->ForcingFunction()(x, f);
//                
//            }
//            
//            
//            Solout[0] = f[0]; // fx
//            Solout[1] = f[1]; // fy
//            Solout[2] = f[2]; // fz
//        }
//            break;
//            
//        case 3: //v_exact
//        {
//            TPZVec<STATE> sol(4,0.0);
//            if(this->HasExactSol()){
//                TPZVec<STATE> x(3,0.);
//                x=datavec[vindex].x;
//                this->fExactSol(x, sol, gradu); // @omar::check it!
//
//            }
//            Solout[0] = sol[0]; // vx
//            Solout[1] = sol[1]; // vy
//            Solout[2] = sol[2]; // vz
//         
//        }
//            break;
//            
//        case 4: //p_exact
//        {
//            TPZVec<STATE> sol(4,0.0);
//            if(this->HasExactSol()){
//                TPZVec<STATE> x(3,0.),xrot(3,0.);
//                x=datavec[pindex].x;
//
//                this->fExactSol(x, sol, gradu); // @omar::check it!
//            }
//            Solout[0] = sol[3]; // px
//            
//        }
//            break;
//
//        case 5: //div
//        {
//            STATE Div=0.;
//            for(int i=0; i<3; i++) {
//                Div+=dsolxy(i,i);
//            }
//            Solout[0] = Div;
//            
//        }
//            break;
//
//        case 6: //norm of tensor
//        {
//            TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
//            TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
//            
//            TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
//            TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
//            
//            int nshapeV;
//            nshapeV = datavec[vindex].fVecShapeIndex.NElements();
//            
//            int normvecRows = datavec[vindex].fDeformedDirections.Rows();
//            int normvecCols = datavec[vindex].fDeformedDirections.Cols();
//            TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
//            TPZManVector<TPZFNMatrix<4,REAL>,18> GradNormalvec(18);
//            
//            STATE asd1 = 0., asd2 = 0.,asd3 = 0., asd4 = 0.;
//            if (datavec[vindex].fNeedsDeformedDirectionsFad) {
//                for (int e = 0; e < normvecRows; e++) {
//                    for (int s = 0; s < normvecCols; s++) {
//                        Normalvec(e,s)=datavec[vindex].fDeformedDirectionsFad(e,s).val();
//                    }
//                }
//                
//                for (int s = 0; s < normvecCols; s++) {
//                    TPZFNMatrix<4,REAL> Grad0(3,3,0.); // 2x2
//                    Grad0(0,0)=datavec[vindex].fDeformedDirectionsFad(0,s).fastAccessDx(0);
//                    Grad0(0,1)=datavec[vindex].fDeformedDirectionsFad(0,s).fastAccessDx(1);
//                    Grad0(1,0)=datavec[vindex].fDeformedDirectionsFad(1,s).fastAccessDx(0);
//                    Grad0(1,1)=datavec[vindex].fDeformedDirectionsFad(1,s).fastAccessDx(1);
//                    GradNormalvec[s] = Grad0;
//                    //Grad0.Print(std::cout);
//                }
//            
//            }else{
//                Normalvec=datavec[vindex].fDeformedDirections;
//            }
//            
//            TPZFMatrix<STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);
//            TPZFNMatrix<4,STATE> GradSol(3,3,0.),GradSolt(3,3,0.),DuSol(3,3,0),GradVi(3,3,0.),GradVit(3,3,0.),Dui(3,3,0.);
//            STATE normDu = 0.;
//            STATE normDuSol = 0.;
//            
//            
//            GradSol = datavec[vindex].dsol[vindex];
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradSolt(e,f) = GradSol(f,e);
//                }
//            }
//            
//            DuSol = GradSolt + GradSol;
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    normDuSol += DuSol(e,f)*DuSol(e,f);
//                }
//            }
//
//            
//            
//            for(int i = 0; i < nshapeV; i++)
//            {
//                int iphi = datavec[vindex].fVecShapeIndex[i].second;
//                int ivec = datavec[vindex].fVecShapeIndex[i].first;
//
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradVi(e,f) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
//                        GradVit(f,e) = Normalvec(e,ivec)*dphiVx(f,iphi)+GradNormalvec[ivec](e,f)*phiV(iphi,0);
//                    }
//                }
//                
//                
//                
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        Dui(e,f)= 0.5 * (GradVi(e,f) + GradVit(e,f));
//                    }
//                }
//                
//               
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        normDu += Dui(e,f)*Dui(e,f);
//                    }
//                }
//
//            }
//
////            GradVi.Print(std::cout);
////            Dui.Print(std::cout);
//            
////            std::cout<<datavec[0].xParametric<<std::endl;
////            std::cout<<datavec[0].x<<std::endl;
////            Normalvec.Print(std::cout);
////            datavec[0].fDeformedDirectionsFad.Print(std::cout);
//            //std::cout<<GradNormalvec<<std::endl;
//            
//            
//            Solout[0] = normDuSol;
//            
//        }
//            break;
//
//        case 7: //PressureCDG is the pressure p (from N-S equation)
//        {
//            Solout[0] = p_h[0]-(v_h[0]*v_h[0]+v_h[1]*v_h[1]+v_h[2]*v_h[2])*0.5;
//        }
//            break;
//        case 8: //p_exact_CDG
//        {
//            TPZVec<STATE> sol(4,0.0);
//            if(this->HasExactSol()){
//                TPZVec<STATE> x(3,0.),xrot(3,0.);
//                x=datavec[pindex].x;
//
//                this->fExactSol(x, sol, gradu); // @omar::check it!
//            }
//            Solout[0] = sol[3]-(sol[0]*sol[0]+sol[1]*sol[1]+sol[2]*sol[2])*0.5; // px
//
//        }
//            break;
//
//        case 9: //Vorticity (multiplicar 1/2???)
//        {
//            Solout[0] = dsolxy(0,1)-dsolxy(1,0);
//        }
//        break;
//            
//        case 10:
//        {
//            TPZFNMatrix<10,STATE> gradUn = datavec[vindex].dsol[0];
//            TPZFNMatrix<9,STATE> DUn_j(3,3,0.), sigma(3,3,0.);
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                     DUn_j(e,f)= 0.5 * (gradUn(e,f) + gradUn(f,e));
//                    sigma(e,f) = 2.*fViscosity*DUn_j(e,f);
//                }
//                sigma(e,e) -= p_h[0];
//            }
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    Solout[e*3+f] = sigma(e,f);
//                }
//            }
//        }
//            break;
//
//        default:
//        {
//            std::cout  << " Var index not implemented " << std::endl;
//            DebugStop();
//        }
//    }
//}
//
//////////////////////////////////////////////////////////////////////
//
//// Divergence on master element
//void TPZNavierStokesMaterial::ComputeDivergenceOnMaster(TPZVec<TPZMaterialDataT<STATE>> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
//
//{
//    int ublock = 0;
//    
//    // Getting test and basis functions
//    TPZFNMatrix<100,REAL> phiuH1        = datavec[ublock].phi;   // For H1  test functions Q
//    TPZFNMatrix<300,REAL> dphiuH1       = datavec[ublock].fDPhi; // Derivative For H1  test functions
//    TPZFNMatrix<300,REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
//    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
//    TPZFNMatrix<9,STATE> graduMaster;
//    gradu.Transpose();
//    
//    TPZFNMatrix<660> GradphiuH1;
//    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
//    
//    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
//    
//    DivergenceofPhi.Resize(nphiuHdiv,1);
//    
//    REAL JacobianDet = datavec[ublock].detjac;
//    
//    TPZFNMatrix<9,REAL> Qaxes = datavec[ublock].axes;
//    TPZFNMatrix<9,REAL> QaxesT;
//    TPZFNMatrix<9,REAL> Jacobian = datavec[ublock].jacobian;
//    TPZFNMatrix<9,REAL> JacobianInverse = datavec[ublock].jacinv;
//    
//    TPZFNMatrix<9,REAL> GradOfX;
//    TPZFNMatrix<9,REAL> GradOfXInverse;
//    TPZFNMatrix<9,REAL> VectorOnMaster;
//    TPZFNMatrix<9,REAL> VectorOnXYZ(3,1,0.0);
//    Qaxes.Transpose(&QaxesT);
//    QaxesT.Multiply(Jacobian, GradOfX);
//    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
//    
//    TPZFMatrix<STATE> GradOfXInverseSTATE(GradOfXInverse.Rows(), GradOfXInverse.Cols());
//    for (unsigned int i = 0; i < GradOfXInverse.Rows(); ++i) {
//        for (unsigned int j = 0; j < GradOfXInverse.Cols(); ++j) {
//            GradOfXInverseSTATE(i,j) = GradOfXInverse(i,j);
//        }
//    }
//    
//    int ivectorindex = 0;
//    int ishapeindex = 0;
//    
//    {
//        for (int iq = 0; iq < nphiuHdiv; iq++)
//        {
//            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
//            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
//            
//            for (int k = 0; k < 3; k++) {
//                VectorOnXYZ(k,0) = datavec[ublock].fDeformedDirections(k,ivectorindex);
//            }
//            
//            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
//            VectorOnMaster *= JacobianDet;
//            
//            /* Contravariant Piola mapping preserves the divergence */
//            for (int k = 0; k < fDimension; k++) {
//                DivergenceofPhi(iq,0) +=  dphiuH1(k,ishapeindex)*VectorOnMaster(k,0);
//            }
//            
//        }
//    }
//    
//    return;
//    
//}
//
//
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::Write(TPZStream &buf, int withclassid) const{
//    
//    TPZMaterial::Write(buf, withclassid);
//    
//    
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::Read(TPZStream &buf, void *context) {
//    
//    TPZMaterial::Read(buf, context);
//    
//}
//
//////////////////////////////////////////////////////////////////////
//
//void TPZNavierStokesMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
//    
//    
//    TPZFMatrix<REAL> &dphiV = dataV.dphix;
//    
//    const int dim = this->Dimension();
//    
//    GradPhi.clear();
//    GradPhi.resize(dim);
//    
//    //for each shape
//    for(int shape = 0; shape < dphiV.Rows(); shape++){
//        
//        TPZFMatrix<REAL> GPhi(dim,dim,0.);
//        
//        for(int i = 0; i < dim; i++){
//            
//            for(int j = 0; j < dim; j++){
//                
//                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
//                
//            }//j
//        }//i
//        
//        GradPhi[shape] = GPhi;
//        
//    }//shape
//    
//}
//
//// Contricucao dos elementos internos
//void TPZNavierStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
//
//    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
//    this->Contribute(datavec, weight, ek_fake, ef);
//
//    return;
//    
//}
//
//// Contricucao dos elementos internos
//void TPZNavierStokesMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
//    
//    // vindex is 0
//    // pindex is 1
//    const int vindex = this->VIndex();
//    const int pindex = this->PIndex();
//
//    // tensorizing the scalar shape functions (in order to use Taylor Hood for instance)
//    // but how does one compute the gradient of these vectors? -> take the gradients of the H1 functions
//    
//    int normvecRows = datavec[vindex].fDeformedDirections.Rows();
//    int normvecCols = datavec[vindex].fDeformedDirections.Cols();
//    
//    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
//    if (datavec[vindex].fVecShapeIndex.size() == 0) {
//        FillVecShapeIndex(datavec[vindex]);
//    }
//    Normalvec=datavec[vindex].fDeformedDirections;
//
//    // substitute by an enumerated variable
//    REAL factorM = 1.;
//    REAL factorMk = 1.;
//    // Setting the phis
//    // V -> we shouldnt be using these variables
//    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
//    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
//    // P
//    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
//    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
//
//
//    int nshapeV, nshapeP;
//    nshapeP = phiP.Rows();
//    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
//
//
//    // a vector of gradients of the H(div) vectors
//    TPZManVector<TPZFNMatrix<9,REAL>,18> GradNormalvec(normvecCols);
//    for (int i=0; i<normvecRows; i++) {
//        GradNormalvec[i].Redim(3,3);
//    }
//
//    if (datavec[vindex].fNeedsDeformedDirectionsFad) {
//        for (int e = 0; e < normvecRows; e++) {
//            for (int s = 0; s < normvecCols; s++) {
//                Normalvec(e,s)=datavec[vindex].fDeformedDirectionsFad(e,s).val();
//            }
//        }
//        //if(fDimension != 2) DebugStop();
//        TPZFNMatrix<9,REAL> Grad0(3,3,0.); // 2x2
//        for (int s = 0; s < normvecCols; s++) {
//            for (int i = 0; i < fDimension; i++) {
//                for (int j = 0; j < fDimension; j++) {
//                    Grad0(i,j)=datavec[vindex].fDeformedDirectionsFad(i,s).fastAccessDx(j);
//                }
//            }
//            GradNormalvec[s] = Grad0;
//        }
//    }
//
//    TPZVec<STATE> Force(3,0.), Force_rot(3,0.);
//    if(f_problemtype==TStokesAnalytic::EBrinkman){
//        Force.resize(4);
//        Force[3] = 0.;
//    }
//
//
//    TPZFMatrix<STATE> phiVi(3,1,0.0),phiVj(3,1,0.0);
//
//    TPZFNMatrix<100,STATE> divphi;
//    TPZFNMatrix<40,STATE> div_on_master;
//    TPZFNMatrix<10,STATE> dsolVec = datavec[vindex].dsol[0];
//    // u_n is the solution at the previous iteration
//
//    TPZManVector<STATE,3> u_n    = datavec[vindex].sol[0];
//    STATE p_n                  = datavec[pindex].sol[0][0];
//    STATE divsol = datavec[vindex].divsol[0][0];
//
//    TPZFNMatrix<10,STATE> dsolp_n = datavec[pindex].dsol[0];
//
//    fDeltaT = f_sim_data->GetTimeStep();
//
//    TPZFNMatrix<10,STATE> gradUn(dsolVec.Rows(),dsolVec.Cols()), grad_axes;
//    grad_axes = datavec[vindex].axes;
//    grad_axes.Resize(dsolVec.Rows(),dsolVec.Cols());
//
//    //std::cout<<grad_axes<<std::endl;
//    gradUn=dsolVec;
//    //TPZAxesTools<REAL>::Axes2XYZ( dsolVec,gradUn, grad_axes);
//
////    if (fSpace==1) {
////        datavec[0].ComputeFunctionDivergence();
////    }
//
//    int64_t global_point_index = datavec[0].intGlobPtIndex;
//
//    if(this->HasForcingFunction()){
//        TPZFMatrix<STATE> gradu;
//        TPZVec<STATE> x(3,0.),xrot(3,0.);
//        x=datavec[vindex].x;
//        this->ForcingFunction()(x, Force);
//    }
//
//    /// LastState indicates that the memory needs to be updated
//    if(fState==ELastState&&fDeltaT>0){
//
//        GetMemory()->operator[](global_point_index).Set_u_last(u_n);
//        return;
//
//    }
///// fDeltaT != 0 means evolution in time. ef needs to accumulate part of the solution at time N
//    if(fState==ECurrentState&&fDeltaT>0){
//
//        // Get the pressure at the integrations points
//        TPZManVector<STATE,3> u_last;
//        TPZNSMemory &lastStep_mem = GetMemory()->operator[](global_point_index);
//        u_last    = lastStep_mem.u_last();
//        // std::cout << u_last << std::endl;
//        // std::cout << u_n<< std::endl;
//
//        for(int ivec = 0; ivec < nshapeV; ivec++ ) {
//            for (int e = 0; e < 3; e++) {
//                phiVi(e, 0) = Normalvec(e, ivec);
//            }
//            STATE phi_dot_Ulast = 0.0; // phi * u_{n-1} / Dt
//            for (int e=0; e<3; e++) {
//                phi_dot_Ulast += phiVi(e)*u_last[e];
//            }
//            ef(ivec) += weight * phi_dot_Ulast /fDeltaT;
//        }
//
//    }
//
//    REAL norm_beta = 0.,norm_un;
//    REAL h_size = datavec[1].HSize;
//    TPZManVector<STATE> beta(3,0.);
//    TPZFNMatrix<9,STATE> grad_beta(3,3,0.), Beta_vec(3,1,0.),un_vec(3,1,0.);
//
//    for (int iv = 0; iv < 3; ++iv) {
//        un_vec(iv,0)=u_n[iv];
//    }
//    norm_un = Norm(un_vec);
//    /// beta contains the known velocity vector, so that the problem becomes linear
//    if (f_problemtype==TStokesAnalytic::EOseenCDG||f_problemtype==TStokesAnalytic::EOseen) {
//        if (HasExactSol()) {
//            TPZVec<STATE> x(3, 0.), xrot(3, 0.);
//            x = datavec[vindex].x;
//            this->fExactSol(x, beta, grad_beta);
//        }
//
//        for (int iv = 0; iv < 3; ++iv) {
//            Beta_vec(iv,0)=beta[iv];
//        }
//        norm_beta = Norm(Beta_vec);
//    }
//
//
//    for(int i = 0; i < nshapeV; i++ )
//    {
////        phiVi : value of the test function
////        GradVi : gradient of the test function
//        // GradVit : transpose of the gradient of the test function (overkill?)
//        int ivec = i;
//        TPZFNMatrix<9,STATE> GradVi(3,3,0.),GradVit(3,3,0.),Dui(3,3,0.);
//        for (int e=0; e<3; e++) {
//            phiVi(e,0) = Normalvec(e,ivec);
//            for (int f=0; f<3; f++) {
//                GradVi(e,f) = GradNormalvec[ivec](e,f);
//                GradVit(f,e) = GradVi(e,f);
//            }
//        }
//        // symmetric gradient of the test function
//        for (int e=0; e<3; e++) {
//            for (int f=0; f<3; f++) {
//                Dui(e,f)= 0.5 * (GradVi(e,f) + GradVi(f,e));
//            }
//        }
//
//        if(fDeltaT>0){
//            STATE phi_dot_Un = 0.0; // - phi * u_{n} / Dt
//            for (int e=0; e<3; e++) {
//                phi_dot_Un += phiVi(e)*u_n[e];
//            }
//            ef(i) += - weight * phi_dot_Un /fDeltaT;
//        }
//
//        /// we are working in residual form (even if the problem is linear)
//        if (f_problemtype==TStokesAnalytic::EBrinkman){
//
//            STATE B_phi_dot_Un = 0.0; // - coefB * phi * u_{n}
//            for (int e=0; e<3; e++) {
//                B_phi_dot_Un += phiVi(e)*u_n[e];
//            }
//            ef(i) += - weight *fcBrinkman * B_phi_dot_Un;
//        }
//
//        STATE divui = 0.;
//        divui = datavec[vindex].divphi(ivec,0);
////        divui = Tr( GradVi ); //datavec[0].divphi(i);
//        //divui = datavec[0].divphi(i);
//
//
//        // phiVi is the vector test function -> should only compute if there is a forcing function
//        STATE phi_dot_f = 0.0, un_dot_phiV = 0.0; // f - Source term
//        for (int e=0; e<3; e++) {
//            phi_dot_f += phiVi(e)*Force[e];
//        }
//        ef(i) += weight * phi_dot_f;
//
//        // computing the residual of the viscous term
//        // gradUn is the gradient of the solution
//        // Dui is the symmetric gradient of the test function
//        // change names? very confusing!
//        STATE A_term_f = 0.; // A - Flux term
//        TPZFNMatrix<9,STATE> DUn_j(3,3,0.);
//        for (int e=0; e<3; e++) {
//            for (int f=0; f<3; f++) {
//                 DUn_j(e,f)= 0.5 * (gradUn(e,f) + gradUn(f,e));
//            }
//        }
//        A_term_f = Inner(Dui, DUn_j);
//
//        ef(i) += 2. * fViscosity * weight * (-A_term_f);
//
//        // why two negatives ? Following formulation signs
//        // p_n is the value of the pressure
//        // divui is the divergence of the test function
//        STATE B_term_f = 0.; // B - Mixed term
//        B_term_f = - p_n * divui;
//        ef(i) += weight * (-B_term_f);
//
//        //std::cout<<ef<<std::endl;
//        TPZFNMatrix<9,STATE> GradUn_phiU(3,1,0.), GradV_phiU(3,1,0.);
//        /// maybe a boolean to indicate a stabilized formulation?
//        STATE C_term_f = 0., Stab_term_f = 0.; // C - Trilinear terms
//        if (f_problemtype==TStokesAnalytic::ENavierStokes) {
//
//            // this is the original Navier Stokes inertia term
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradUn_phiU(e,0) += gradUn(e,f)*u_n[f];
//                }
//            }
//
//            C_term_f = InnerVec(GradUn_phiU, phiVi);
//
//            ef(i) += -weight * C_term_f;
//
//
//            //Stabilization-NS:
////            for (int e=0; e<3; e++) {
////                for (int f=0; f<3; f++) {
////                    GradV_phiU(e,0) += GradVi(e,f)*u_n[f];
////                }
////            }
////            Stab_term_f = InnerVec(GradUn_phiU, GradV_phiU);
////            ef(i) += -0.5*h_size*(1./norm_un)*weight * Stab_term_f; //nsok
//
//        }
//
//        TPZFNMatrix<9,STATE> GradUnTr_phiU(3,1,0.),GradVTr_phiU(3,1,0.),WUn_Un(3,1,0.),WV_Un(3,1,0.),WV_V(3,1,0.);
//        TPZFNMatrix<9,STATE> GradV_phiV(3,1,0.),GradVTr_phiV(3,1,0.);
//        // multiplying by the anti-symmetric of the gradient of u_n
//        if (f_problemtype==TStokesAnalytic::ENavierStokesCDG) {
//
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradUn_phiU(e,0) += gradUn(e,f)*u_n[f];
//                    GradUnTr_phiU(e,0) += gradUn(f,e)*u_n[f];
//                }
//            }
//
//            C_term_f = InnerVec(GradUn_phiU, phiVi) - InnerVec(GradUnTr_phiU, phiVi);
//
//            ef(i) += -weight * C_term_f;
//
//        }
//
//
//        TPZFNMatrix<9,STATE> GradUn_Beta(3,1,0.),GradV_Beta(3,1,0.),GradVTr_Beta(3,1,0.),GradUnTr_Beta(3,1,0.);
//        TPZFNMatrix<9,STATE> GradBeta_Beta(3,1,0.);
//        if (f_problemtype==TStokesAnalytic::EOseen) {
//            // this is an adjustment for a Oseen version of Navier Stokes?
//            // beta is the value of the velocity of the exact solution
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradUn_Beta(e,0) += gradUn(e,f)*beta[f];
//                }
//            }
//
//            C_term_f = InnerVec(GradUn_Beta, phiVi);
//            ef(i) += -weight * C_term_f;
//
//            //Stabilization-Oseen:
////            for (int e=0; e<3; e++) {
////                for (int f=0; f<3; f++) {
////                    GradV_phiU(e,0) += GradVi(e,f)*u_n[f];
////                    GradV_Beta(e,0) += GradVi(e,f)*beta[f];
////                    GradBeta_Beta(e,0)+= grad_beta(e,f)*beta[f];
////                }
////            }
////            STATE Stab_term_f = InnerVec(GradBeta_Beta, GradV_Beta);
////            ef(i) += -0.5*h_size*(1./norm_beta)*weight * Stab_term_f;
//        }
//
//        TPZFNMatrix<9,STATE> WUn_Beta(3,1,0.),WV_Beta(3,1,0.),WBeta_Beta(3,1,0.);
//        TPZFNMatrix<9,STATE> GradTrBeta_Beta(3,1,0.);
//
//
//        if (f_problemtype==TStokesAnalytic::EOseenCDG) {
//
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradUn_Beta(e,0) += gradUn(e,f)*beta[f];
//                    GradUnTr_Beta(e,0) += gradUn(f,e)*beta[f];
//                }
//            }
//
//            C_term_f = InnerVec(GradUn_Beta, phiVi) - InnerVec(GradUnTr_Beta, phiVi);
//            ef(i) += -weight * C_term_f;
//
//        }
//
//
//        // A, C e D - velocity X velocity
//        for(int j = 0; j < nshapeV; j++){
//            int jvec = j;
//
//            for (int e=0; e<3; e++) {
//                phiVj(e,0) = Normalvec(e,jvec);
//            }
//
//            if(fDeltaT>0){
//                STATE Transient_term = InnerVec(phiVi, phiVj);
//                ek(i,j) += weight * Transient_term / fDeltaT;  //phiV * phiU / Dt
//            }
//
//            if (f_problemtype==TStokesAnalytic::EBrinkman){
//                STATE Brinkman_term = fcBrinkman * InnerVec(phiVi, phiVj); // - coefB * phiU * phiV
//                ek(i,j) += weight  * Brinkman_term;
//            }
//
//            TPZFNMatrix<9,STATE> GradVj(3,3,0.),GradVjt(3,3,0.),Duj(3,3,0.);
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    GradVj(e,f) = GradNormalvec[jvec](e,f);
//                    GradVjt(f,e) = GradVj(e,f);
//                }
//            }
//            for (int e=0; e<3; e++) {
//                for (int f=0; f<3; f++) {
//                    Duj(e,f)= 0.5 * (GradVj(e,f) + GradVj(f,e));
//                }
//            }
//            STATE A_term = Inner(Dui, Duj);
//            STATE viscc = fViscosity;
//            ek(i,j) += 2. * weight * fViscosity * A_term;  // A - Bilinear gradV * gradU
//
//
//            if (f_problemtype==TStokesAnalytic::ENavierStokes) {
//
//                TPZFNMatrix<9,STATE> GradU_Un(3,1,0.), GradUTr_Un(3,1,0.);
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradU_Un(e,0) += GradVj(e,f)*u_n[f]; //Oseen eqs
//                    }
//                }
//
//                TPZFNMatrix<9,STATE> GradUn_phiVj(3,1,0.);
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradUn_phiVj(e,0) += gradUn(e,f)*phiVj(f,0);
//                    }
//                }
//
//                STATE C_term = InnerVec(GradU_Un, phiVi);
//
//                STATE C_term_2 = InnerVec(GradUn_phiVj, phiVi);
//
//                // THIS IS WRONG!! EITHER TRANSPOSE OR NOT
//                ek(i,j) += weight * C_term;  // C - Trilinear terms
//
//                ek(i,j) += weight * C_term_2;  // C - Trilinear terms
//
//
//                //Stabilization-NS-CDG:
////                TPZFNMatrix<9,STATE> GradVi_phiVj(3,1,0.);
////
////                for (int e=0; e<3; e++) {
////                    for (int f=0; f<3; f++) {
////                        GradVi_phiVj(e,0) += GradVi(e,f)*phiVj(f,0);
////                    }
////                }
////
////                STATE Stab_term = InnerVec(GradU_Un, GradV_phiU)+InnerVec(GradUn_phiVj, GradV_phiU)+InnerVec(GradUn_phiU, GradV_phiU);
////                ek(i,j) += 0.5*h_size*(1./norm_un)*weight * Stab_term; //nsok
//
//            }
//
//            if (f_problemtype==TStokesAnalytic::ENavierStokesCDG) {
//
//                TPZFNMatrix<9,STATE> GradU_Un(3,1,0.), GradUTr_Un(3,1,0.);
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradU_Un(e,0) += GradVj(e,f)*u_n[f]; //Oseen eqs Obs:GradU = GradVj = Grad Du
//                        GradUTr_Un(e,0)+= GradVj(f,e)*u_n[f];
//                    }
//                }
//
//                TPZFNMatrix<9,STATE> GradUn_phiVj(3,1,0.),GradUnTr_phiVj(3,1,0.);
//
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradUn_phiVj(e,0) += gradUn(e,f)*phiVj(f,0);
//                        GradUnTr_phiVj(e,0) += gradUn(f,e)*phiVj(f,0);
//                    }
//                }
//
//                STATE C_term = InnerVec(GradU_Un, phiVi) - InnerVec(GradUTr_Un, phiVi);
//
//                STATE C_term_2 = InnerVec(GradUn_phiVj, phiVi) - InnerVec(GradUnTr_phiVj, phiVi);
//
//                ek(i,j) += weight * C_term;  // C - Trilinear terms
//
//                ek(i,j) += weight * C_term_2;  // C - Trilinear terms
//
//
//
//            }
//
//
//            if (f_problemtype==TStokesAnalytic::EOseen) {
//
//                TPZFNMatrix<9,STATE> GradU_Beta(3,1,0.);
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradU_Beta(e,0) += GradVj(e,f)*beta[f]; //Oseen eqs
//                    }
//                }
//
//                STATE C_term = InnerVec(GradU_Beta, phiVi);
//
//                ek(i,j) += weight * C_term;  // C - Trilinear terms
//
//                //Stabilization-Oseen:
////                STATE Stab_term = InnerVec(GradU_Beta, GradV_Beta);
////                ek(i,j) += 0.5*h_size*(1./norm_beta)*weight * Stab_term;
//
//            }
//
//            if (f_problemtype==TStokesAnalytic::EOseenCDG) {
//
//                TPZFNMatrix<9,STATE> GradU_Beta(3,1,0.), GradUTr_Beta(3,1,0.);
//                for (int e=0; e<3; e++) {
//                    for (int f=0; f<3; f++) {
//                        GradU_Beta(e,0) += GradVj(e,f)*beta[f]; //Oseen eqs
//                        GradUTr_Beta(e,0) += GradVj(f,e)*beta[f];
//                    }
//                }
//
//                STATE C_term = InnerVec(GradU_Beta, phiVi) - InnerVec(GradUTr_Beta, phiVi);
//
//                ek(i,j) += weight * C_term;  // C - Trilinear terms
//
//                //Stabilization-OseenCDG:
////                TPZFNMatrix<9,STATE> WU_Beta(3,1,0.);
////                for (int e=0; e<3; e++) {
////                    WU_Beta(e,0) = (1./2.)*(GradU_Beta(e,0)-GradUTr_Beta(e,0));
////                }
////                STATE Stab_term = 2.*InnerVec(WU_Beta, WV_Beta);
////                ek(i,j) += 0.5*h_size*(1./norm_beta)*weight * Stab_term; //ekocdg
////
////                for (int j = 0; j < nshapeP; j++) {
////
////                    TPZFNMatrix<9,STATE> GradPj(3,1,0.);
////                    for (int e=0; e<3; e++) {
////                        GradPj(e,0) = phiV(e,j);
////                    }
////
////                    STATE PStab_term = 0.;
////                    PStab_term = 0.5*h_size*(1./norm_beta)*weight *InnerVec(phiVj, WV_Beta);
////                    // Matrix B
////                    ek(i, nshapeV+j) += PStab_term; //zxzxzxzx
////                    // Matrix B^T
////                    ek(nshapeV+j,i) += PStab_term;  //zxzxzxzx
////                }
////
//            }
//
//
//        }
//
//        // B - pressure and velocity
//        for (int j = 0; j < nshapeP; j++) {
//
//
//            STATE B_term = 0.;
//            B_term = (-1.) * weight * phiP(j,0) * divui;
//            // Matrix B
//            ek(i, nshapeV+j) += B_term; // B - Bilinear div v * p
//            // Matrix B^T
//            ek(nshapeV+j,i) += B_term;  // Bt - Bilinear div u * q
//
//        }
//
//    }
//
//    // VERIFY!!!
//    for (int i = 0; i < nshapeP; i++) {
//
//        STATE B_term_f = 0.; // B - Mixed term
//        B_term_f = - phiP(i,0)*divsol;
//        ef(i+nshapeV) += weight * (-B_term_f);
//
//        if (f_problemtype==TStokesAnalytic::EBrinkman) {
//            STATE Brinkman_source = 0.; // B - Mixed term
//            Brinkman_source = phiP(i,0)*Force[3]; //gsource
//            ef(i+nshapeV) += -weight*Brinkman_source;
//        }
//
//    }
//
//
//    // Preparação para formulação MHM :
//    if (datavec.size()>2) {
//        
//        TPZFMatrix<REAL> &phigM = datavec[2].phi;
//        TPZFMatrix<REAL> &phipM = datavec[3].phi;
//        STATE g_average = datavec[2].sol[0][0];
//        STATE p_average = datavec[3].sol[0][0];
//        
//        // matrix D - pressure and average-pressure
//        for (int j = 0; j < nshapeP; j++) {
//            
//            STATE fact = (1.) * weight * phiP(j,0) * phigM(0,0);
//            // Matrix D
//            ek(nshapeV+nshapeP, nshapeV+j) += fact;
//            // Matrix D^T
//            ek(nshapeV+j,nshapeV+nshapeP) += fact;
//            ef(nshapeV+j,0) -= weight*phiP(j,0) *g_average;
//        }
//        
//        ef(nshapeV+nshapeP) -= (p_n-p_average) * weight;
//        ef(nshapeV+nshapeP+1) -= -g_average * weight;
//        
//        // matrix E - injection and average-pressure
//
//        STATE factG = (1.) * weight * phigM(0,0) * phipM(0,0);
//        // Matrix E
//        ek(nshapeV+nshapeP+1, nshapeV+nshapeP) += -factG;
//        // Matrix E^T
//        ek(nshapeV+nshapeP,nshapeV+nshapeP+1) += -factG;
//    }
//
//    if (datavec.size()>4) {
//        
//        TPZFMatrix<REAL> &phigM0 = datavec[4].phi;
//        TPZFMatrix<REAL> &phipM0 = datavec[5].phi;
//        STATE g_average = datavec[4].sol[0][0];
//        STATE p_average = datavec[5].sol[0][0];
//
//        // matrix D0 - pressure and average-pressure
//        for (int j = 0; j < nshapeP; j++) {
//            
//            STATE fact0 = (1.) * weight * phiP(j,0) * phigM0(0,0);
//            // Matrix D
//            ek(nshapeV+nshapeP+2, nshapeV+j) += fact0;
//            // Matrix D^T
//            ek(nshapeV+j,nshapeV+nshapeP+2) += fact0;
//            
//        }
//        //DebugStop();
//        // matrix E0 - injection and average-pressure
//        ef(nshapeV+nshapeP+2) -= (p_n-p_average) * weight;
//        ef(nshapeV+nshapeP+3) -= -g_average * weight;
//
//        STATE factG0 = (1.) * weight * phigM0(0,0) * phipM0(0,0);
//        // Matrix E
//        ek(nshapeV+nshapeP+3, nshapeV+nshapeP+2) += -factG0;
//        // Matrix E^T
//        ek(nshapeV+nshapeP+2,nshapeV+nshapeP+3) += -factG0;
//        
//    }
//
//
//#ifdef PZDEBUG
//    if(0)
//    {
//        std::ofstream fileEK("FileEKContribute.txt");
//        ek.Print("stiff = ",fileEK,EMathematicaInput);
//        
//        std::ofstream fileEF("FileEFContribute.txt");
//        ef.Print("rhs = ",fileEF,EMathematicaInput);
//    }
//#endif
//    
//}
//
//
//void TPZNavierStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
//    
//    STATE rhsnorm = Norm(ef);
//    if(isnan(rhsnorm))
//    {
//        std::cout << "ef  has norm " << rhsnorm << std::endl;
//    }
//
//#ifdef PZDEBUG
//    //2 = 1 Vel space + 1 Press space
//    int nref =  datavec.size();
//    if (nref != 2 ) {
//        std::cout << " Erro. The size of the datavec is different from 2 \n";
//        DebugStop();
//    }
//#endif
//    
//    
//    const int vindex = this->VIndex();
//    const int pindex = this->PIndex();
//    
//    if (datavec[vindex].fVecShapeIndex.size() == 0) {
//        // I don't understand this, so DebugStop()...
//        DebugStop();
//        FillVecShapeIndex(datavec[vindex]);
//    }
//  
//    // Setting the phis
//    // V
//    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
//    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
//    // P
//    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
//    
//    // Getting the linear combination or finite element approximations
//    
//    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
//    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
//    
//    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
////    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
//    
//    dphiVx = dphiV;
//    
//    int nshapeV, nshapeP;
//    nshapeP = phiP.Rows();
// //   nshapeV = phiV.Rows()*NStateVariables();
//    nshapeV = datavec[vindex].fDeformedDirections.Cols();//datavec[vindex].fVecShapeIndex.NElements();
//    
//    int normvecRows = datavec[vindex].fDeformedDirections.Rows();
//    int normvecCols = datavec[vindex].fDeformedDirections.Cols();
//    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
//    
////    if (datavec[vindex].fNeedsDeformedDirectionsFad==false) {
//        Normalvec=datavec[vindex].fDeformedDirections;
////    }else{
////        for (int e = 0; e < normvecRows; e++) {
////            for (int s = 0; s < normvecCols; s++) {
////                Normalvec(e,s)=datavec[vindex].fDeformedDirectionsFad(e,s).val();
////            }
////        }
////    }
//    
////    Normalvec.Print(std::cout);
//    
//    if (fSpace==1) {
//        nshapeV = nshapeV/2.;
//    }
//
//
//    int gy=v_h.size();
//    
//    TPZFNMatrix<9,STATE> phiVi(fDimension,1,0.),phiVni(1,1,0.), phiVj(fDimension,1,0.),phiVnj(1,1,0.), phiPi(fDimension,1),phiPj(fDimension,1);
//    
//    TPZManVector<STATE,3> v_2=bc.Val2();
//    TPZFNMatrix<3,STATE> v_1=bc.Val1();
//    STATE p_D = bc.Val1()(0,0);
//    
//    switch (bc.Type()) {
//        case 0: //Dirichlet for continuous formulation
//        {
//            
//            TPZFMatrix<STATE> gradu(3,3,0.);
//            TPZManVector<STATE> vbc(4,0.);
//            TPZFMatrix<STATE> Du(3,3,0.),Dun(3,1,0.);
//
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> vbc(4,0.);
//                TPZFMatrix<STATE> gradu(3,3,0.);
//                bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
//                v_2[0] = vbc[0];
//                v_2[1] = vbc[1];
//                v_2[2] = vbc[2];
//                p_D = vbc[3];
//            }
//            
//            if(fSpace==1){
//
//                
//                for(int i = 0; i < nshapeV; i++ )
//                {
//
//                    //Adaptação para Hdiv
//
//                    TPZManVector<REAL> n = datavec[0].normal;
//
//                    REAL vh_n = v_h[0];
//                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
//
//                    ef(i,0) += -weight * fBigNumber * (vh_n - v_n) * phiV(i,0);
//
//                    for(int j = 0; j < nshapeV; j++){
//
//                        ek(i,j) += weight * fBigNumber * phiV(j,0) * phiV(i,0);
//                        
//                    }
//                    
//
//
//                }
//                
//                
//                
//            }else{
//                
//                for(int i = 0; i < nshapeV; i++ )
//                {
//                    int ivec = i;
//                    
//                    for (int e=0; e<fDimension; e++) {
//                        phiVi(e,0)=Normalvec(e,ivec);
//                    }
//                    
//                    
//                    //Adaptação para Hdiv
//                    
//                    STATE factef=0.0;
//                    for(int is=0; is<gy ; is++){
//                        factef += -1.0*(v_h[is] - v_2[is]) * phiVi(is,0);
//                    }
//                    
//                    ef(i,0) += weight * fBigNumber * factef;
//                    
//                    for(int j = 0; j < nshapeV; j++){
//                        int jvec = j;
//                        
//                        for (int e=0; e<fDimension; e++) {
//                            phiVj(e,0)=Normalvec(e,jvec);
//                        }
//                        
//                        //Adaptação para Hdiv
//                        
//                        STATE factek = 0.0;
//                        for(int is=0; is<gy ; is++){
//                            factek += phiVj(is,0) * phiVi(is,0);
//                        }
//                        
//                        ek(i,j) += weight * fBigNumber * factek;
//                        
//                    }
//                    
//                }
//            }
//            
//        }
//            break;
//            
//        case 1: //Neumann for continuous formulation
//        {
//            
//            
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> vbc(4,0.);
//                TPZFMatrix<STATE> gradu(3,3,0.);
//                bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
//                v_2[0] = vbc[0];
//                v_2[1] = vbc[1];
//                v_2[2] = vbc[2];
//                p_D = vbc[3];
//            }
//            
//            
//            for(int i = 0; i < nshapeV; i++ )
//            {
//                int ivec = i;
//                
//                for (int e=0; e<fDimension; e++) {
//                    phiVi(e,0)=Normalvec(e,ivec);
//                }
//                
//                TPZManVector<REAL> n = datavec[vindex].normal;
//                
//                TPZFNMatrix<9,STATE> pn(fDimension,1);
//                
//                
//                for (int f=0; f<fDimension; f++) {
//                    pn(f,0)=n[f]*v_1(0,0);
//                }
//                
//                //Adaptação para Hdiv
//                
//                STATE factef=0.0;
//                for(int is=0; is<gy ; is++){
//                    factef += (pn(is,0))* phiVi(is,0);
//                }
//                
//                ef(i,0) += weight * factef;
//                
//            }
//            
//            
//        }
//            
//            
//            
//            break;
//            
//            
//        case 2: //Condição Penetração
//        {
//            
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> vbc(4,0.);
//                TPZFMatrix<STATE> gradu(3,3,0.);
//                bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
//                v_2[0] = vbc[0];
//                v_2[1] = vbc[1];
//                v_2[2] = vbc[2];
//                p_D = vbc[3];
//            }
//
//            
//            if(fSpace==1){
//
//                
//                for(int i = 0; i < nshapeV; i++ )
//                {
//                    
//                    //Adaptação para Hdiv
//                    
//                    TPZManVector<REAL> n = datavec[0].normal;
//                    
//                    REAL vh_n = v_h[0];
//                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
//                    
//                    ef(i,0) += -weight * fBigNumber * (vh_n - v_n) * phiV(i,0);
//                    
//                    for(int j = 0; j < nshapeV; j++){
//                        
//                        ek(i,j) += weight * fBigNumber * phiV(j,0) * phiV(i,0);
//                        
//                    }
//                    
//                }
//                
//                
//            }else{
//                
//                
//                
//                
//                TPZManVector<REAL> n = datavec[0].normal;
//                TPZManVector<REAL> t(2);
//                t[0]=-n[1];  //oioioio
//                t[1]=n[0];
//                
//                
//                
//                //Componente normal -> imposta fortemente:
//                
//                for(int i = 0; i < nshapeV; i++ )
//                {
//                    
//                    int ivec = i;
//                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
//                    
//                    
//                    for (int e=0; e<fDimension; e++) {
//                        phiVi(e,0)=Normalvec(e,ivec);
//                        phiVni(0,0)+=phiVi(e,0)*n[e];
//                        phiVti(0,0)+=phiVi(e,0)*t[e];
//                    }
//                    
//                    REAL vh_n = v_h[0];
//                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
//                    
//                    ef(i,0) += -weight * fBigNumber * (vh_n-v_n) * (phiVni(0,0));
//                    
//                    
//                    for(int j = 0; j < nshapeV; j++){
//                        
//                        int jvec = j;
//                        
//                        TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
//                        
//                        for (int e=0; e<fDimension; e++) {
//                            phiVj(e,0)=Normalvec(e,jvec);
//                            phiVnj(0,0)+=phiVj(e,0)*n[e];
//                            phiVtj(0,0)+=phiVj(e,0)*t[e];
//                            
//                        }
//                        
//                        ek(i,j) += weight * fBigNumber * phiVni(0,0) * phiVnj(0,0);
//                        
//                    }
//                    
//                }
//                
//                
//                
//            }
//            
//        }
//            break;
//            
//            
//        case 3: //Contribuicao ponto no x
//        {
//            
//            REAL p_D = v_2[0];
//            
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> pbc(1);
//                TPZFNMatrix<4,STATE> val;
//                bc.ForcingFunctionBC()(datavec[vindex].x,pbc,val);
//                p_D = pbc[0];
//                
//            }
//            
//            TPZManVector<REAL> n = datavec[0].normal;
//            
//            REAL phiVi_n;
//            
//            
//            // THIS IS VERY STRANGE!! WHY ITERATE ONLY FROM ONE TO THREE?
//            // changed - corrected?
//            for(int i = 0; i < nshapeV; i++ )
//            {
//                phiVi_n = 0.0;
//                
//                int ivec = i;
//                
//                for (int e=0; e<fDimension; e++) {
//                    phiVi(e,0)=Normalvec(e,ivec);
//                    phiVi_n += phiVi(e,0)*n[e];
//                }
//                
//                ef(i*2,0) += -1.0*weight * p_D * phiVi_n;
//                
//                
//            }
//            
//            
//        }
//            break;
//            
//        case 4: //Contribuicao ponto no y
//        {
//            
//            REAL p_D = v_2[0];
//            
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> pbc(1);
//                TPZFNMatrix<4,STATE> val;
//                bc.ForcingFunctionBC()(datavec[vindex].x,pbc,val);
//                p_D = pbc[0];
//                
//            }
//            
//            TPZManVector<REAL> n = datavec[0].normal;
//            
//            REAL phiVi_n;
//            
//            
//            for(int i = 0; i < nshapeV; i++ )
//            {
//                phiVi_n = 0.0;
//                
//                int ivec = i;
//                
//                for (int e=0; e<fDimension; e++) {
//                    phiVi(e,0)=Normalvec(e,ivec);
//                    phiVi_n += phiVi(e,0)*n[e];
//                }
//                
//                ef(i*2+1,0) += -1.0*weight * p_D * phiVi_n;
//                
//                
//            }
//            
//            
//        }
//            break;
//            
//        case 5: //Ponto pressao
//        {
//           
//            //return;
//            p_D = bc.Val2()[0];
//            
//            
//            for(int i = 0; i < nshapeP; i++ )
//            {
//                
//                
//                ef(i) += 1.0 * p_D * phiP(i,0);
//                
//                for(int j = 0; j < nshapeP; j++){
//                    
//                    ek(i,j) += 1.0 * (phiP(i,0) * phiP(j,0));
//                    
//                }
//                
//            }
//            
//        }
//            break;
//            
//        case 6: //Pressao Dirichlet
//        {
//            
//            if(bc.HasForcingFunctionBC())
//            {
//                TPZManVector<STATE> vbc(3);
//                TPZFNMatrix<4,STATE> val;
//                bc.ForcingFunctionBC()(datavec[pindex].x,vbc,val);
//                v_2[0] = vbc[0];
//                v_2[1] = vbc[1];
//                p_D  = vbc[2]*0.;
//                
//                
//            }
//            
//            //pressao
//            
//            for(int i = 0; i < nshapeP; i++ )
//            {
//                
//                
//                ef(i) += -weight * fBigNumber * (p_h[0] - p_D) * phiP(i,0);
//                
//                for(int j = 0; j < nshapeP; j++){
//                    
//                    ek(i,j) += weight * fBigNumber * (phiP(i,0) * phiP(j,0));
//                    
//                }
//                
//            }
//            
//        }
//            
//            break;
//            
//            
//        default:
//        {
//            std::cout << "Boundary not implemented " << std::endl;
//            DebugStop();
//        }
//            break;
//    }
//    
//    
//    if(isnan(rhsnorm))
//    {
//        std::cout << "ef  has norm " << rhsnorm << std::endl;
//    }
//
//    
//#ifdef PZDEBUG
//    if(0)
//    {
//        std::ofstream fileEK("FileEKContributeBC.txt");
//        std::ofstream fileEF("FileEFContributeBC.txt");
//        ek.Print("stiff = ",fileEK,EMathematicaInput);
//        ef.Print("force = ",fileEF,EMathematicaInput);
//    }
//#endif
//
//}
//
//////////////////////////////////////////////////////////////////////
//template <typename TVar>
//TVar TPZNavierStokesMaterial::Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
//    
//    //inner product of two tensors
//    
//    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
//        DebugStop();
//    }
//    
//    
//#ifdef DEBUG
//    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
//        DebugStop();
//    }
//#endif
//    
//    TVar Val = 0;
//    
//    for(int i = 0; i < S.Cols(); i++){
//        for(int j = 0; j < S.Cols(); j++){
//            Val += S(i,j)*T(i,j);
//        }
//    }
//    
//    return Val;
//    
//}
//
//
//////////////////////////////////////////////////////////////////////
//STATE TPZNavierStokesMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
//    
//    //inner product of two vectors
//    
//    
//#ifdef DEBUG
//    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
//        DebugStop();
//    }
//#endif
//    
//    STATE Val = 0;
//    
//    for(int j = 0; j < S.Cols(); j++){
//        for(int i = 0; i < S.Rows(); i++){
//            Val += S(i,j)*T(i,j);
//        }
//    }
//    
//    return Val;
//    
//}
//
//
//
//////////////////////////////////////////////////////////////////////
//
//STATE TPZNavierStokesMaterial::Tr( TPZFMatrix<REAL> &GradU ){
//    
//#ifdef DEBUG
//    if( GradU.Rows() != GradU.Cols() ) {
//        DebugStop();
//    }
//#endif
//    
//    STATE Val = 0.;
//    
//    for(int i = 0; i < GradU.Rows(); i++){
//        Val += GradU(i,i);
//    }
//    
//    return Val;
//}
//
//
///// transform a H1 data structure to a vector data structure
//void TPZNavierStokesMaterial::FillVecShapeIndex(TPZMaterialData &data)
//{
//    int64_t nshape = data.phi.Rows();
//    data.fDeformedDirections.Redim(fDimension,nshape*fDimension);
//    for (int d=0; d<fDimension; d++) {
//        for (int i=0; i<data.phi.Rows(); i++) {
//            data.fDeformedDirections(d,i*fDimension+d) = data.phi(i,0);
//        }
//    }
//    TPZFNMatrix<60,REAL> dphix(3,nshape*fDimension);
//    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, dphix, data.axes);
//    DebugStop();
//}
//
//
//
//void TPZNavierStokesMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
//{
//    
//    TPZManVector<STATE,4> sol_exact(3);
//    TPZFNMatrix<16> dsol_exact(3,3);
//    if(!HasExactSol()) DebugStop();
//    fExactSol(data[0].x,sol_exact,dsol_exact);
//    errors.Resize(NEvalErrors());
//    errors.Fill(0.0);
//    TPZManVector<STATE> Velocity(3,0.), Pressure(3,0.);
//
//    this->Solution(data,VariableIndex("V"), Velocity);
//    this->Solution(data,VariableIndex("P"), Pressure);
//
//    int vindex = this->VIndex();
//    int pindex = this->PIndex();
//    
//    TPZFMatrix<REAL> dudx(3,3);
//    TPZFMatrix<STATE> &dsol = data[vindex].dsol[0];
//    TPZFMatrix<STATE> &dsolp = data[pindex].dsol[0];
//    //std::cout<<dsol<<std::endl;
//    
//    //Adaptação feita para Hdiv
//    dsol.Resize(3,3);
//    
//    TPZFNMatrix<2,STATE> dsolxy(3,3,0.), dsolxyp(3,1,0.);
//    dsolxy = dsol;
//    if (fSpace!=1) {
//        TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data[vindex].axes);
//    }
//   // TPZAxesTools<STATE>::Axes2XYZ(dsolp, dsolxyp, data[pindex].axes);
//    
//    dsolxyp = dsolp;
//
////    std::cout<<Velocity<<std::endl;
////    std::cout<<sol_exact<<std::endl;
//
//    int shift = 3;
//    // velocity - erro norma L2
//    STATE diff, diffp;
//    errors[0] = 0.;
//    for(int i=0; i<3; i++) {
//        diff = Velocity[i] - sol_exact[i];
//        errors[0]  += diff*diff;
//    }
//
//    // velocity - erro divergence
//
//    STATE Div_exact=0., Div=0.;
//    for(int i=0; i<3; i++) {
//        Div_exact+=dsol_exact(i,i);
//        Div+=dsolxy(i,i);
//    }
//
//    diff = Div-Div_exact;
//    errors[1]  = diff*diff;
//
//    // pressure - eror em norma L2
//    diffp = Pressure[0] - sol_exact[3];
//    errors[2]  = diffp*diffp;
//
//
//    //For couplings:
//    if(f_problemtype==TStokesAnalytic::EBrinkman){
//        for (int i = 0; i < 3; ++i) {
//            errors[3+i]=errors[i];
//            errors[i]=0.;
//        }
//    }
//
//    if(f_problemtype==TStokesAnalytic::ENavierStokesCDG||f_problemtype==TStokesAnalytic::EOseenCDG){
//        STATE diffp_CDG = 0.;
//        diffp = (Pressure[0]-(Velocity[0]*Velocity[0]+Velocity[1]*Velocity[1]+Velocity[2]*Velocity[2])*0.5) - (sol_exact[3]-(sol_exact[0]*sol_exact[0]+sol_exact[1]*sol_exact[1]+sol_exact[2]*sol_exact[2])*0.5);
//        errors[shift+2]  = diffp*diffp;
//    }
//
//    
//}
//
