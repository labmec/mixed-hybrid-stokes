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

#include "ProblemData.h"
#include "SimpleExample.h"
#include "TPZMeshOperator.h"
#include <pzskylstrmatrix.h>

int main()
{
    TPZFMatrix<std::complex<double>> C(3,3,1.0);
    TPZFMatrix<std::complex<double>> A(2,1,1.0);
    TPZFMatrix<std::complex<double>> B(2,2,2.0);

    C.AddContribution(1, 1, A, true, B, false, 1.0);

    C.Print("C", std::cout, EMathematicaInput);

    return 0;
}