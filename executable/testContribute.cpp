
#include <iostream>
#include <functional>
#include <cmath>
#include <TPZGenGrid3D.h>
#include <TPZGenGrid2D.h>
#include <DarcyFlow/TPZDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZLinearAnalysis.h>
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZMultiphysicsCompMesh.h"
#include <pzmultiphysicscompel.h>
#include <TPZCompElH1.h>
#include <pzelchdiv.h>
#include <pzinterpolationspace.h>
#include "TPZSimpleTimer.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZVTKGenerator.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "pzskylstrmatrix.h"
#include "TPZVTKGeoMesh.h"
#include <TPZMaterialDataT.h>
#include "TPZNullMaterial.h"
#include "TPZNavierStokesMaterial.h"
#include "pzgeoquad.h"
#include <pzshapequad.h>
#include <TPZAnalyticSolution.h>
#include "TPZElementMatrixT.h"

using namespace std;

enum EnumMatids {
  EMatId = 1,  // Material ID for the domain
  EBottom = 2,
  ERight = 3,
  ETop = 4,
  ELeft = 5
};

TPZGeoMesh* createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX) {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId);  // 1 is the material ID for the domain
  generator.SetBC(gmesh, 4, EBottom);   // Set boundary conditions
  generator.SetBC(gmesh, 5, ERight);
  generator.SetBC(gmesh, 6, ETop);
  generator.SetBC(gmesh, 7, ELeft);

  return gmesh;
}

TPZGeoMesh* createRectangularGmesh() {
  
  TPZGeoMesh *gmesh = new TPZGeoMesh;

  TPZManVector<REAL,3> coord(3,0.);

  TPZGeoNode nod0(0, coord, *gmesh);
  coord[0] = 1.;
  TPZGeoNode nod1(1, coord, *gmesh);
  coord[1] = 1.;
  TPZGeoNode nod2(2, coord, *gmesh);
  coord[0] = 0.;
  TPZGeoNode nod3(3, coord, *gmesh);
  coord[0] = 2.;
  coord[1] = 0.;
  TPZGeoNode nod4(4, coord, *gmesh);
  coord[1] = 1.;
  TPZGeoNode nod5(5, coord, *gmesh);
  gmesh->NodeVec().Resize(6);
  gmesh->NodeVec()[0] = nod0;
  gmesh->NodeVec()[1] = nod1;
  gmesh->NodeVec()[2] = nod2;
  gmesh->NodeVec()[3] = nod3;
  gmesh->NodeVec()[4] = nod4;
  gmesh->NodeVec()[5] = nod5;

  // --- First quadrilateral element ---
  TPZManVector<int64_t,4> nodeIndexes(4);
  nodeIndexes[0] = 0;
  nodeIndexes[1] = 1;
  nodeIndexes[2] = 2;
  nodeIndexes[3] = 3;

  int64_t index = 0;
  TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // --- Second quadrilateral element ---
  nodeIndexes[0] = 1;
  nodeIndexes[1] = 4;
  nodeIndexes[2] = 5;
  nodeIndexes[3] = 2;
  index = 1;
  TPZGeoEl *gel2 = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // Creating BCs
  int iddiri = 2;
  // gel->CreateBCGeoEl(4, EBottom); 
  // gel->CreateBCGeoEl(5, ERight);
  // gel->CreateBCGeoEl(6, ETop);
  // gel->CreateBCGeoEl(7, ELeft);
  gel->CreateBCGeoEl(4,EBottom);
  gel2->CreateBCGeoEl(4,EBottom);
  gel2->CreateBCGeoEl(5,ERight);
  gel2->CreateBCGeoEl(6,ETop);
  gel->CreateBCGeoEl(6,ETop);
  gel->CreateBCGeoEl(7,ELeft);

  gmesh->SetDimension(2);  
  gmesh->BuildConnectivity();

  return gmesh;
}

TPZCompMesh* createCompMeshFlux(TPZGeoMesh *gmesh) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->SetAllCreateFunctionsHDiv(); // Set all element creation functions to HDiv

  // Add materials (weak formulation)
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());  
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions
  int diritype = 0, neumanntype = 1, mixedtype = 2;
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  // Create boundary conditions
  // val2[0] = 0.;
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // val2[0] = 1.;
  bcond = mat->CreateBC(mat, ERight, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  // val2[0] = 0.;
  bcond = mat->CreateBC(mat, ETop, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  // val2[0] = 3.;
  bcond = mat->CreateBC(mat, ELeft, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
}

// -------------------------------------------- PRESSURE COMP MESH -----------------------------------------------------
TPZCompMesh* createCompMeshPressure(TPZGeoMesh *gmesh) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  const int pOrder = 1;
  cmesh->SetDefaultOrder(pOrder);
  if (pOrder < 1) {
    cmesh->SetAllCreateFunctionsDiscontinuous();
  } else {
    cmesh->SetAllCreateFunctionsContinuous(); // Set all element creation functions to Continuous H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
  }

  // Add materials (weak formulation)
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());
  cmesh->InsertMaterialObject(mat);

  // Set up the computational mesh
  cmesh->AutoBuild();

  //! ASK THIS
  int ncon = cmesh->NConnects();
  const int lagLevel = 1; // Lagrange multiplier level
  for(int i=0; i<ncon; i++)
  {
      TPZConnect &newnod = cmesh->ConnectVec()[i]; 
      newnod.SetLagrangeMultiplier(lagLevel);
  }

  return cmesh;
}

//-------------------------------------------- MULTIPHYSICS COMP MESH ----------------------------------------------------
TPZMultiphysicsCompMesh* createMPCompMesh(TPZGeoMesh *gmesh, TPZCompMesh* cmeshFlux, TPZCompMesh* cmeshPressure) {
  TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

 // Add materials (weak formulation)
  const REAL visc =  1.;
  //TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
  TPZNavierStokesMaterial *mat = new TPZNavierStokesMaterial(EMatId, gmesh->Dimension(), visc);
  // TPZStokesMaterial *mat = new TPZStokesMaterial(EMatId, gmesh->Dimension(), visc);
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions
  int diritype = 0, neumanntype = 1, mixedtype = 2;
  TPZManVector<REAL,1> val2(1,3.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix

  // Create boundary conditions
  val2[0] = 0.;
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumanntype, val1, val2);


  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = mat->CreateBC(mat, ERight, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  val2[0] = 0.;
  bcond = mat->CreateBC(mat, ETop, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  val2[0] = 3.;
  bcond = mat->CreateBC(mat, ELeft, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Add the computational meshes to the multiphysics mesh
  TPZManVector<TPZCompMesh *,2> cmeshes(2);
  cmeshes[0] = cmeshFlux; // Flux mesh
  cmeshes[1] = cmeshPressure; // Pressure mesh

  TPZManVector<int> active(cmeshes.size(),1);    
  cmesh->BuildMultiphysicsSpace(active, cmeshes);

  return cmesh;
}

//--------------------------------------- TESTING CONTRIBUTE NAVIER-STOKES-OSEEN -------------------------------------------
void BetaOseen(const TPZVec<REAL> &x, TPZVec<REAL> &beta) {
  // Define the beta function for the Oseen problem
  beta[0] = x[0]*x[1];  // Example: constant value for beta
  beta[1] = x[0]+x[1];
  if (beta.size() > 2){
    beta[3] = 1.;
  }
}

void TestContributeNS(TPZGeoMesh *gmesh, TPZCompMesh* cmeshFlux, TPZCompMesh* cmeshPressure, TPZMultiphysicsCompMesh* cmesh, TPZElementMatrixT<STATE>& ek, TPZElementMatrixT<STATE>& ef){ //TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef

  //TODO Preliminary
  int dim = cmesh->Dimension(); // Problem dimension
  const REAL visc =  1.; // Setting viscosity 
  TPZNavierStokesMaterial *NSmat = dynamic_cast<TPZNavierStokesMaterial*>(cmesh->FindMaterial(EMatId)); // Initializing new material from the computational mesh material
  int64_t NCompEl = cmesh->NElements(); // Number of computational elements in comp mesh
  std::cout << "Number of computational elements in comp mesh: " << NCompEl << std::endl;

  //TODO Getting multiphysics element
  TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>* CompElMulti = dynamic_cast<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>*>(cmesh->ElementVec()[1]); // Getting (multiphysics) elem [0] from comp mesh
  /* dynamic_cast<new_type>(expression): converts data from one type to another type at runtime. Mainly used in inherited class hierarchies for safely 
  casting the base (polymorphic, with at least one virtual function) class pointer or reference to derived class 
  (called downcasting). Example:
  Base* basePtr = new Derived();
  Derived* derivedPtr = dynamic_cast<Derived*>(basePtr)   (downcasting);*/
  int64_t NCompMeshes = CompElMulti->NMeshes(); // Number of computational meshes associated with the multiphysics elem (flux and pressure in this case)
  std::cout << "Number of meshes: " << NCompMeshes << std::endl;

  // CompEl<pzgeom::TPZGeoQuad>* CompElFlux = dynamic_cast<CompEl<pzgeom::TPZGeoQuad>*>(cmeshFlux->ElementVec()[0]);
  // CompEl<pzgeom::TPZGeoQuad>* CompElPressure = dynamic_cast<CompEl<pzgeom::TPZGeoQuad>*>(cmeshPressure->ElementVec()[0]);

  // TPZInterpolationSpace* CompElFlux = dynamic_cast<TPZInterpolationSpace*>(cmeshFlux->ElementVec()[0]);
  // TPZInterpolationSpace* CompElPressure = dynamic_cast<TPZInterpolationSpace*>(cmeshPressure->ElementVec()[0]);

  // TPZCompElHDiv<pzshape::TPZShapeQuad>* CompElFluxHDiv = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad>*>(cmeshFlux->ElementVec()[0]); 
  // TPZCompElH1<pzshape::TPZShapeQuad>* CompElPressureH1 = dynamic_cast<TPZCompElH1<pzshape::TPZShapeQuad>*>(cmeshPressure->ElementVec()[0]);

  //TODO Configuring Material Data (like in TPZMultiphysicsCompEl<TGeometry>::CalcStiffT)
  TPZManVector<TPZMaterialDataT<STATE>,5> datavec(NCompMeshes); // Initializing datavec with number of meshes associated
  TPZVec<REAL> qsi(2); // Coordinates of one integration point in master element (randomly chosen)
  qsi[0] = 0.7;
  qsi[1] = 0.3; 
  REAL weight = 1.; // Weight of integration rule (randomly chosen)

  // CompElFluxHDiv->InitMaterialData(datavec[0]);
  // CompElPressureH1->InitMaterialData(datavec[1]);
  // CompElFluxHDiv->ComputeRequiredData(datavec[0],qsi);
  // CompElPressureH1->ComputeRequiredData(datavec[1],qsi);

  CompElMulti->InitMaterialData(datavec); // Initialize a material data vector and its attributes, TPZMultiphysicsCompEl<TGeometry>::InitMaterialDataT
  TPZVec<TPZTransform<> > trvec;
  CompElMulti->AffineTransform(trvec); // Map a parameter point in the multiphysic elem to a parameter point in the master element, TPZMultiphysicsCompEl<TGeometry>::AffineTransform
  CompElMulti->ComputeRequiredData(qsi, trvec, datavec); // Compute and fill data with requested attributes, TPZMultiphysicsElement::ComputeRequiredDataT
  CompElMulti->InitializeElementMatrix(ek,ef); // Initialize element matrix in which is computed CalcStiff, TPZMultiphysicsCompEl<TGeometry>::InitializeElementMatrix

  // CompElMulti->Print(std::cout);
  datavec[0].Print(std::cout);
  datavec[1].Print(std::cout);
  
  //TODO Setting an Analytical Solution and Forcing Function 
  // To be attributed to Beta Function
  TStokesAnalytic *analyticSol = new TStokesAnalytic(); // Getting an exact solution already implemented 
  analyticSol->fvisco = visc;
  analyticSol->fvelocity = 1.;
  analyticSol->fconstPressure = 1.;
  analyticSol->fProblemType = TStokesAnalytic::EOseen; // Used in Pablo's Doctor //{EStokes, ENavierStokes, EOseen, ENavierStokesCDG, EOseenCDG, EBrinkman}
  analyticSol->fExactSol = TStokesAnalytic::EKovasznay; // Used in Pablo's Doctor //{ENone, ECavity,  EKovasznay, EKovasznayCDG, ESinCos, ENoFlow, ESinCos3D, EPconst, EObstacles, EOneCurve , ESinCosBDS, ESinCosBDS3D, EGatica3D, ECouplingSD, ECouplingNSD, EVugs2D, EVugs3D, EInfiltrationNS, EConstantFlow, ECouetteFlow, EPoisFlow, ETaylorCouette, EElbow, EPaperComp, ESimpleF}
  analyticSol->fDimension = dim;

  NSmat->SetExactSol(analyticSol->ExactSolution(), 2);
  // NSmat->SetForcingFunction(analyticSol->ForceFunc(), 2); //! DO IT

  //TODO Setting Beta Function 
  NSmat->setBetaFunc(BetaOseen);
  // analyticSol->Solution(datavec[0].x, Beta, gradBeta); // Setting Beta as an analytical solution

  // TODO Contribute
  NSmat->Contribute(datavec, weight, ek.fMat, ef.fMat);

}

int main (int argc, char * const argv[]) {

  const int nthreads = 0;

  // ------- Create geometric mesh -------
  const bool isUseGenGrid = false; // Use TPZGenGrid2D for mesh generation
  TPZGeoMesh* gmesh = nullptr;
  if (isUseGenGrid) {
    const int neldiv = 1;
    gmesh = createMeshWithGenGrid({neldiv, neldiv}, {-1., -1.}, {1., 1.});
  } else {
    gmesh = createRectangularGmesh();
  }
  // gmesh->Print(std::cout);
  // std::ofstream out("geomesh.vtk");
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // ------- Create computational mesh -------
  // Create flux Hdiv cmesh
  TPZCompMesh* cmeshFlux = createCompMeshFlux(gmesh);
  // Create pressure L2 cmesh
  TPZCompMesh* cmeshPressure = createCompMeshPressure(gmesh);
  // Create multiphysics cmesh

  TPZMultiphysicsCompMesh* cmesh = createMPCompMesh(gmesh, cmeshFlux, cmeshPressure);  
  // cmesh->Print(std::cout);


  // Initializing element Stiffness matrix and Load vector
  // TPZFMatrix<STATE> ke;
  // TPZFMatrix<STATE> fe;
  TPZElementMatrixT<STATE> ke;
  TPZElementMatrixT<STATE> fe;
  // ke.Reset();
  // fe.Reset();

  TestContributeNS(gmesh, cmeshFlux, cmeshPressure, cmesh, ke, fe);
  std::stringstream sout;
  ke.fMat.Print("kel", std::cout, EMathematicaInput);
  // ke.Print(std::cout);



//   // ------- Create Analysis object -------
//   TPZLinearAnalysis an(cmesh);
//   // TPZSSpStructMatrix<STATE> matsp(cmesh);
//   TPZSkylineStructMatrix<STATE> matsp(cmesh);
//   matsp.SetNumThreads(nthreads); // Set number of threads for assembly
//   an.SetStructuralMatrix(matsp);
//   TPZStepSolver<STATE> step;
//   step.SetDirect(ELDLt); // Set direct solver type
//   an.SetSolver(step);
//   an.Run();
//   an.Solution().Print("Solution");

//   // Post-process the results
//   const std::string plotfile = "postproc";  // sem o .vtk no final
//   constexpr int vtkRes{0};

//   TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
//   auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

//   vtk.Do();

//   delete cmesh;
//   delete gmesh;

  return 0;
}
