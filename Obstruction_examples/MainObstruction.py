"""
This script is a new version of how to create pipe modules with obstructions
"""
#%% ****************** 
#   IMPORTED MODULES
#   ******************
import gmsh

#%% ****************** 
#   IMPORTED CLASSES
#   ******************
from TPZMeshModeling import TPZMeshModeling
from TPZSimpleObstruction import TPZSimpleObstruction
from TPZCrossObstruction import TPZCrossObstruction
from TPZRandomObstruction import TPZRandomObstruction
from TPZMultipleObstruction import TPZMultipleObstruction
from TPZSemiArcObstruction import TPZSemiArcObstruction
from TPZNoObstruction import TPZNoObstruction

#%% ****************** 
#     JSON DATA
#   ******************
file_name = "Reference"

mm = 1e-3
cm = mm*10

# "Input for module generation"
length = 40*cm
radius = 45*mm

lc = 1e-2

# "Input for obstruction generation"
obstruction = 80*mm/2

json_data = {
        "MeshName": "../Meshes/"+file_name,
        "CreateMsh": False,
        "HdivType": 1,
        "VelpOrder": 2,
        "TracpOrder": 1,
        "Dim": 3,
        "Resolution": 0,
        "StaticCondensation": True,
        "isAxisymmetric": 0,
        "HasAnalyticSolution": False,
        "Domain": [
            {
                "name": "Domain",
                "matID": 1,
                "viscosity": 1
            }
        ],
        "NormalBoundary": [
            {
                "name": "VelIn",
                "type": 2,
                "value": -1,
                "matID": 2
            },
            {
                "name": "PressOut",
                "type": 2,
                "value": 0,
                "matID": 3
            },
            {
                "name": "NoPenetration",
                "type": 0,
                "value": 0,
                "matID": 5
            }
        ],
        "TangentialBoundary": [
            {
                "name": "NoSlip",
                "type": 1,
                "value": [0,0,0],
                "matID": 4
            }
        ],
        "AxisymmetryDomain": [
        ],
        "AxisymmetryBoundary":[
        ],
        "LambdaID": 10,
        "InterfaceID": 20,
        "AxiLambdaID": 30,
        "AxiInterfaceID": 40,
        "FluxInterfaceID": 50,
        "ObstructionID": 100
    }
#%% ****************** 
#     MAIN FUNCTION
#   ******************
def main()->None:
    """
    Main function
    """
    TPZMeshModeling.Begin()

    TPZMeshModeling.TurnOnLabels('surface', 'volume')
    TPZMeshModeling.TurnOnRepresentation('surfaces')
    
    TPZMeshModeling.TurnOnNormals()
    # TPZMeshModeling.TurnOnTangents()

    # "Input for mesh generation"
    mesh_dim = 3
    circle = ('Circular', {'radius': radius})

    # "Creating the obstructions"
    circular = True

    modules = []
    if circular:
        # modules.append(TPZSimpleObstruction(length, lc, circle, obstruction))
        # modules.append(TPZCrossObstruction(length, lc, circle, 2e-2, 2e-2, 1e-3))
        # modules.append(TPZMultipleObstruction(length, lc, circle, 5e-3, 3e-2))
        modules.append(TPZNoObstruction(_length = length, _lc = lc, _module_typology = circle))
        modules.append(TPZNoObstruction(_length = length, _lc = lc, _module_typology = circle))
    
    "Moving them to the right place"
    module: TPZSimpleObstruction
    for i, module in enumerate(modules):
        module.Move(0, 0, length*i)
    
    # "Joining them all in one"
    gmsh.model.occ.removeAllDuplicates()

    physical_group = [
        [(3, [i + 1 for i, _ in enumerate(modules) ]), 1, "Domain"],
        [(2, [1]), 2, "PressIn"],
        [(2, [6]), 3, "PressOut"],
        [(2, [2, 3, 4, 5, 8, 9, 10, 11, 1 ,6]), 4, "NoSlip"],
        [(2, [2, 3, 4, 5, 8, 9, 10, 11]), 5, "NoPenetration"],
        # [(2, [6]), 100, "Obstruction"]
    ]

    # "Creating the elements"
    TPZMeshModeling.Synchronize()

    TPZMeshModeling.CreatePhysicalGroup(physical_group)

    TPZMeshModeling.CreateMesh(mesh_dim)
    
    TPZMeshModeling.ShowModel()

    TPZMeshModeling.WriteMeshFiles(file_name, ".msh")

    TPZMeshModeling.PrintJson(json_data, file_name)

    TPZMeshModeling.End()

if __name__ == '__main__':
    main()