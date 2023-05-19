#include <iostream>
#include <fstream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>
#include <TPZSSpStructMatrix.h>
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include <TPZVTKGenerator.h>
#include <TPZSimpleTimer.h>

#include"ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("Stokes.cfg");
#endif
    
    
    bool printdata = true;

    std::string filepath = "../DataInput/";
    std::string filename =  "LidDrivenFlow";

    ProblemData simData;
    simData.ReadJson(filepath+filename+".json");

    TPZGeoMesh* gmesh = TPZMeshOperator::CreateGMesh(&simData);
    TPZCompMesh* cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);
    TPZCompMesh* cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);

    TPZMultiphysicsCompMesh* cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh);

    TPZLinearAnalysis an(cmesh_m, true);
    TPZSSpStructMatrix<> strmat(cmesh_m);

    strmat.SetNumThreads(0);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);

    an.Assemble();

    an.Solve();
    
    if(printdata){
        simData.Print();
        TPZMeshOperator::PrintMesh(gmesh, cmesh_v, cmesh_p, cmesh_m);
    }

    //vtk export
    {
        TPZSimpleTimer timer("PostProcess", true);
        
        TPZVTKGenerator vtk(cmesh_m, {"Pressure", "Velocity", "Tension"}, filename,simData.Resolution());
        vtk.Do();
    }
    
    std::cout << "\n\nSimulation finished without errors :) \n\n";
    
	return 0;
}
