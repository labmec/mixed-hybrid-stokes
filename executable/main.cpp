#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzfstrmatrix.h>
#include <pzlog.h>
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include <TPZVTKGenerator.h>
#include <TPZSimpleTimer.h>

#include "ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"
#include <pzskylstrmatrix.h>
#include <pzfstrmatrix.h>
#include <pzstrmatrixot.h>
#include <TPZSpStructMatrix.h>
#include "TPZIdentitySolver.h"
#include <TPZPardisoSolver.h>
#include "TPZGeoMeshTools.h"
#include "pzvec_extras.h"

const int global_nthread = 32;

struct PostProcElData {
  TPZGeoEl* gel;
  TPZManVector<REAL,1> qsi = {-10}; // a value that is not valid
  TPZManVector<REAL,3> x = {0.,0.,0.};
};

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh, std::string filename, ProblemData *problem_data, TPZGeoMesh *gmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name);
void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh, ProblemData *problem_data);
void ComputeIntegralFlux(TPZCompMesh *cmesh, int matid, std::ofstream &out);
void FindElementsAndPtsToPostProc(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, const int matid, const int npts, const REAL z0, const REAL zf);
void PostProcDataForArticle(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, std::ofstream &out);

int main()
{
#ifdef PZ_LOG
//    TPZLogger::InitializePZLOG("log4cxx.cfg");
    TPZLogger::InitializePZLOG();
#endif
    
    bool printdata = false;
    
    std::string filepath = "/home/giavancini/dev/obstructor-analyses/mixed-hybrid-stokes/";
    std::string filename = "Obstructor-10mm";

    ProblemData simData;
    simData.ReadJson(filepath + filename + ".json");
    
    // creating a analytical solution for ExactSol problem
    TStokesAnalytic *flow = new TStokesAnalytic();
    flow->fvisco = simData.DomainVec()[0].viscosity;
    flow->fvelocity = 1.;
    flow->fconstPressure = 1.;
    flow->fExactSol = TStokesAnalytic::ENone;
    flow->fDimension = simData.Dim();
    
    // creating the meshes
    // geometric mesh
    TPZGeoMesh *gmesh = nullptr;
    gmesh = TPZMeshOperator::CreateGMesh(&simData, false);
    
    // velocity computational mesh
    TPZCompMesh *cmesh_v = nullptr;
    cmesh_v = TPZMeshOperator::CreateCMeshV(&simData, gmesh);

    // pressure computational mesh
    TPZCompMesh *cmesh_p;
    cmesh_p = TPZMeshOperator::CreateCmeshP(&simData, gmesh);
    
    // atomic meshes for static condensation
    TPZCompMesh *cmesh_Mp = nullptr; // mean pressure
    TPZCompMesh *cmesh_Mv = nullptr; // and mean flux
    
    if(simData.CondensedElements() && simData.HdivType() != ProblemData::EConstant)
    {
        cmesh_Mv = TPZMeshOperator::CreateCmeshG(&simData, gmesh);
        cmesh_Mp = TPZMeshOperator::CreateCmeshPm(&simData, gmesh);
    }

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = nullptr;
    cmesh_m = TPZMeshOperator::CreateMultiPhysicsMesh(&simData, gmesh, flow);
    
    // static condensation
    if (simData.CondensedElements())
        TPZMeshOperator::CondenseElements(&simData, cmesh_m, gmesh);
        
    TPZLinearAnalysis an(cmesh_m, RenumType::EMetis);
    
    // starting the simulation
    SolveProblem(an, cmesh_m, filename, &simData, gmesh);

    std::ofstream out(filename + "-data.txt");
    ComputeIntegralFlux(cmesh_m, 3, out);
    TPZStack<PostProcElData> postProcData;
    FindElementsAndPtsToPostProc(gmesh, postProcData, simData.DomainVec()[0].matID, 1001, 0.0, 1.0);
    PostProcDataForArticle(gmesh, postProcData, out);
        
    // printing meshes
    if (printdata)
    {
        cmesh_m->ComputeNodElCon();
        TPZMeshOperator::PrintCompMesh(simData.MeshVector());
        TPZMeshOperator::PrintCompMesh(cmesh_m);
        
        TPZMeshOperator::PrintGeoMesh(gmesh);
    }
    
    // Post Process
    PrintResults(an, cmesh_m, &simData, filename);
    
    // Calculating error
    if (flow->fExactSol)
        EvaluateErrors(filename, an, flow, cmesh_m, &simData);
    
    // Deleting pointers
    if (cmesh_m)
        delete cmesh_m;
    
    if (cmesh_v)
        delete cmesh_v;
    
    if (cmesh_p)
        delete cmesh_p;
    
    if (cmesh_Mp)
        delete cmesh_Mp;
    
    if (cmesh_Mv)
        delete cmesh_Mv;
    
    if (gmesh)
        delete gmesh;
    
    std::cout << "\n\nSimulation finished without errors :) \n\n";
    
	return 0;
}

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, std::string filename, ProblemData *problem_data, TPZGeoMesh *gmesh)
{
    TPZSimpleTimer timer;
    
    TPZSSpStructMatrix<> strmat(cmesh_m);
    strmat.SetNumThreads(global_nthread);
        
    //    Applying filters to null normal and tangential velocities
    TPZEquationFilter filter(cmesh_m->NEquations());
    std::set<int64_t> removeEquations;

    //     on the boundary
    TPZMeshOperator::ConfigureBoundaryFilter(gmesh, cmesh_m, problem_data, removeEquations);

    //     between the domain and the obstruction
    if (problem_data->ObstructionID() != - 1)
        TPZMeshOperator::ConfigureObstructionFilter(gmesh, cmesh_m, problem_data, removeEquations);

    filter.ExcludeEquations(removeEquations);
    strmat.EquationFilter() = filter;
    
    an.SetStructuralMatrix(strmat);
    
    // setting solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    // ... the assemble
    std::cout << "Starting assemble...\n";
    TPZSimpleTimer ass_time("Assemble time");
    
    an.Assemble();
    
    std::cout << "Time to assemble = " << ass_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished assemble...\n";
    
    // ... and the solver
    std::cout << "Starting solver...\n";
    std::cout << "Number of Equations " << cmesh_m->NEquations() << std::endl;
    
    TPZSimpleTimer solve_time("Solve time");
    
    an.Solve();
    
    // saving the simulation data (time of execution and number of equations)
    std::cout << "Time to solve = " << solve_time.ReturnTimeDouble()/1000. << " s" << std::endl;
    std::cout << "Finished solver...\n";
    
    std::cout << "Simulation Time: " << timer.ReturnTimeDouble()/1000. << std::endl;
    
    return;
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh_m, ProblemData *problem_data, std::string file_name)
{
    std::cout << "Starting PostProcess" << std::endl;
    TPZSimpleTimer postProc("Post proc time");
    
    // vtk export
    TPZVec<std::string> fields;
    
    if (problem_data->HasAnalyticSolution())
        fields =
        {
            "Pressure",
            "ExactPressure",
//            "ErrorPressure",
            
            "Velocity",
            "ExactVelocity",
//            "ErrorVelocity",
            
            "Stress",
            "ExactStress",
//            "ErrorStress"
        };
    else
        fields = {"Pressure", "Velocity", "Stress"};
    
    TPZVTKGenerator vtk(cmesh_m, fields, file_name, problem_data->Resolution());
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
    
    return;
}

void EvaluateErrors(std::string file_name, TPZLinearAnalysis &an, TPZAnalyticSolution *flow, TPZCompMesh *cmesh_m, ProblemData *problem_data)
{
    std::cout << "Calculating error..." << std::endl;
    
    TPZVec<std::string> fields =
    {
        "p error",
        "p_ex error",
        
        "u error",
        "u_ex error",
        
        "div error",
        "div_ex error",
        
        "sigma error",
        "sigma_ex error",
        
        "dev_sigma error",
        "dev_sigma_ex error"
    };
    
    an.SetExact(flow->ExactSolution());
    an.SetThreadsForError(global_nthread);
    
    TPZMaterial *mat = cmesh_m->FindMaterial(problem_data->DomainVec()[0].matID);
    TPZMatErrorCombinedSpaces<STATE> *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    TPZManVector<REAL, 10> Errors(matError->NEvalErrors());
    
    bool store_errors = false;
    std::ofstream ErrorOut(file_name  + "_Errors.txt");
    
    TPZSimpleTimer post_time("PostProcess time");
    
    an.PostProcessError(Errors, store_errors);
    
    std::cout << "Time PostProc Error = " << post_time.ReturnTimeDouble()/1000 << "s" << std::endl;
    
    ErrorOut << "###### Computed Errors ######" << std::endl;
    
    for (int i = 0; i < Errors.size(); i++)
        ErrorOut << "L2 " << fields[i] << " = " << Errors[i] << std::endl;
}

void ComputeIntegralFlux(TPZCompMesh *cmesh, int matid, std::ofstream& out)
{
    //Getting the equations number of the DOFs associated with the material id
    int64_t nel = cmesh->NElements();
    std::set<int64_t> equations;
    REAL area = 0.;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) continue; //condensed elements and element group have no reference
        if (gel->MaterialId() != matid) continue;
        area += gel->Volume();
        int ncon = cel->NConnects();
        for (int ic = 0; ic < ncon; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int seqnum = c.SequenceNumber();
            int pos = cmesh->Block().Position(seqnum);
            int blocksize = cmesh->Block().Size(seqnum);
            for (int ieq = 0; ieq < blocksize; ieq++) {
                equations.insert(pos+ieq);
            }
        }
    }

    //Sum the flux coefficients
    REAL integralflux = 0.0;
    for (auto eq : equations) {
        TPZFMatrix<STATE>& sol = cmesh->Solution();
        integralflux += sol(eq,0);
    }
    std::cout << "Integral flux = " << integralflux << std::endl;
    out << "Integral flux = " << integralflux << std::endl;
    std::cout << "Area = " << area << std::endl;
    out << "Area = " << area << std::endl;
}

void FindElementsAndPtsToPostProc(TPZGeoMesh *gmesh, TPZStack<PostProcElData> &postProcData, const int matid, const int npts, const REAL z0, const REAL zf)
{
    const REAL dz = (zf - z0) / (npts - 1);

    int64_t InitialElIndex = -1;
    for (auto gel : gmesh->ElementVec())
    {
        if (gel && gel->MaterialId() == matid)
        {
            InitialElIndex = gel->Index();
            break;
        }
    }

    TPZManVector<REAL, 3> xvec(3, 0.);
    xvec[0] = xvec[1] = 0;
    for (int i = 0; i < npts; i++)
    {
        const REAL z = z0 + i * dz;
        xvec[2] = z;
        TPZManVector<REAL, 3> qsi(3, 0.0);
        TPZGeoEl *gel = TPZGeoMeshTools::FindElementByMatId(gmesh, xvec, qsi, InitialElIndex, {matid});
        if (!gel)
        {
            std::cout << "Element not found for z = " << z << std::endl;
            DebugStop();
        }
#ifdef PZDEBUG
        TPZManVector<REAL, 3> xcheck(3);
        gel->X(qsi, xcheck);
        REAL distance = dist(xvec, xcheck);
        if (distance > 1.e-8)
        {
            DebugStop(); // check if the element found is the correct one
        }
#endif
        postProcData.Push({gel, qsi, xvec}); // Creates a struct entry with gel and qsi.
    }
}

void PostProcDataForArticle(TPZGeoMesh *gmesh, TPZStack<PostProcElData> &postProcData, std::ofstream &out)
{
    std::ofstream out_graph("plot_over_line.txt");
    out_graph << std::setprecision(15);
    REAL deltaP = 0;
    for (auto &data : postProcData)
    {
        TPZVec<REAL> qsi(3, 0.0);
        TPZCompEl *cel = data.gel->Reference();
        if (!cel)
            DebugStop();
        TPZMaterial *mat = cel->Material();
        TPZStokesMaterial *stokesmat = dynamic_cast<TPZStokesMaterial *>(mat);
        if (!stokesmat)
            DebugStop();
        const int uind = stokesmat->VariableIndex("Velocity");
        const int pind = stokesmat->VariableIndex("Pressure");
        TPZManVector<STATE, 9> output(3);
        cel->Solution(data.qsi, uind, output);
        const REAL velz = output[2];
        cel->Solution(data.qsi, pind, output);
        const REAL pressure = output[0];

        if (fabs(data.x[2]-0.475) <= 1.e-8)
        {
            deltaP += pressure;
        }
        else if (fabs(data.x[2]-0.525) <= 1.e-8)
        {
            deltaP -= pressure;
        }

        out_graph << data.x[2] << " " << velz << " " << pressure << std::endl;
    }
    std::cout << "DeltaP = " << deltaP << std::endl;
    out << "DeltaP = " << deltaP << std::endl;
}