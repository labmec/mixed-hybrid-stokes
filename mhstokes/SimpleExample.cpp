#include <iostream>
#include "TPZGenGrid2D.h"
#include <TPZGenGrid3D.h>
#include <pzvec.h>
#include <pzgmesh.h>
#include <TPZVTKGeoMesh.h>
#include <TPZVTKGenerator.h>
#include <Poisson/TPZMatPoisson.h>
#include <pzcmesh.h>
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>

void ExampleH1(){
    std::string geofile = "QuadGmesh.vtk";
    std::string compfile = "QuadCmesh.vtk";
    
    constexpr int dim{2};
    
    // The approximation space -> [-1,-1,-1]x[1,1,1]
    // lower corner
    const TPZManVector<REAL, 3> X0 = {0.,0.,0.};
    //upper corner
    const TPZManVector<REAL, 3> X1 = {1.,1.,1.};
    
    // number of divisions in x,y and z directions
    int ndiv = 20;
    
    const TPZManVector<int, 3> division = {ndiv, ndiv, ndiv};
    
    // and finally creates the geometric mesh
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    
    // generates the grid
    if (dim == 2){
        TPZGenGrid2D generator(division, X0, X1);
        
        // Set the 2D grid
        generator.SetElementType(MMeshType::EQuadrilateral);
        
        generator.Read(gmesh);
        
        // Defines the material ID of each corner
        generator.SetBC(gmesh, 4, 2);
        generator.SetBC(gmesh, 5, 2);
        generator.SetBC(gmesh, 7, 2);
        
        generator.SetBC(gmesh, 6, 3);
    } else if (dim == 3){
        TPZGenGrid3D generator(X0, X1, division, MMeshType::ETriangular);
        
        gmesh = generator.BuildVolumetricElements(1);
        gmesh = generator.BuildBoundaryElements(2, 2, 2, 2, 2, 3);
    }
    
    // prints the mesh in VTK format
    std::ofstream outfile(geofile);
    
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile);
    
    // creates the computational mesh from the geometric one
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    
    cmesh -> SetAllCreateFunctionsContinuous();
    
    // Defines the material and boundary conditions
    // first for the element domain
    
    auto *mat = new TPZMatPoisson<STATE>(1, dim);
    cmesh -> InsertMaterialObject(mat);
    
    // and now for the domain boundary
    TPZFMatrix<STATE> val1(1.,1.,0.); // goes in the matrix
    TPZManVector<STATE, 1> val2_e{0}; // goes in the rhs;
    TPZManVector<STATE, 1> val2_t{1}; // different due to the unit velocity vector applied on the top
    
    constexpr int boundType{0};
    
    TPZBndCond* bnd2 = mat-> CreateBC(mat, 2, boundType, val1, val2_e);
    cmesh->InsertMaterialObject(bnd2);
    
    TPZBndCond* bnd3 = mat-> CreateBC(mat, 3, boundType, val1, val2_t);
    cmesh->InsertMaterialObject(bnd3);
    
    // creates the computational mesh
    constexpr int nOrder{2};
    cmesh->SetDefaultOrder(nOrder);
    cmesh->AutoBuild();
    
    std::ofstream coutfile(compfile);
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, coutfile);
    
    // creates an Analysis Object
    TPZLinearAnalysis analysis(cmesh);
    
    //Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    {
      TPZSimpleTimer total("Total");
      {
        TPZSimpleTimer assemble("Assemble");
        //assembles the system
        analysis.Assemble();
      }
      {
        TPZSimpleTimer solve("Solve");
        ///solves the system
        analysis.Solve();
      }
    }
    
    //vtk export
      TPZVec<std::string> scalarVars(1), vectorVars(1);
      scalarVars[0] = "Solution";
      vectorVars[0] = "Derivative";
      analysis.DefineGraphMesh(dim,scalarVars,vectorVars,"poissonSolution.vtk");
      constexpr int resolution{0};
      analysis.PostProcess(resolution);
    
}
