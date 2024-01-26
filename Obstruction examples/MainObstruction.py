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
from TPZNoObstruction import TPZNoObstruction
from TPZCrossObstruction import TPZCrossObstruction
from TPZMultipleObstruction import TPZMultipleObstruction
from TPZSemiArcObstruction import TPZSemiArcObstruction
from TPZRandomObstruction import TPZRandomObstruction

#%% ****************** 
#     JSON DATA
#   ******************

json_data = {
        "MeshName": "../Meshes/SimpleObstruction",
        "CreateMsh": False,
        "HdivType": 1,
        "VelpOrder": 1,
        "TracpOrder": 0,
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
                "name": "PressIn",
                "type": 2,
                "value": -10,
                "matID": 2
            },
            {
                "name": "PressOut",
                "type": 2,
                "value": 10,
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
                "name": "NoTanVel",
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
    # File name
    file_name = "SimpleObstruction"

    gmsh.initialize()

    TPZMeshModeling.TurnOnLabels('surfaces', 'volumes')
    # TPZMeshModeling.TurnOnNumbering('points', 'volumes', 'surfaces')
    
    TPZMeshModeling.TurnOnNormals()
    TPZMeshModeling.TurnOnTangents()

    # "Input for mesh generation"
    lc = 1e-1
    mesh_dim = 3

    # "Input for module generation"
    length = .5
    radius = .5
    circle = ('Circular', {'radius': radius})
    rec = ('Rectangular', {'dx': 2*radius, 'dy':2*radius})
 
    # "Input for obstruction generation"
    r = .1
    width = .3
    height = .3
    obstruction = .25

    # "Creating the obstructions"
    circular = True

    modules = []
    if circular:
        # modules.append(TPZCrossObstruction(_length = length ,_lc = lc, _module_typology = circle, _radius = r/2, _obstruction_width = width, _obstruction_height = height))
        modules.append(TPZMultipleObstruction(length, lc, circle, obstruction/3, .3))
        modules.append(TPZSimpleObstruction(length, lc, circle, obstruction))
        # modules.append(TPZSemiArcObstruction(length, lc, circle, obstruction))
        # modules.append(TPZRandomObstruction(length, lc, circle, obstruction/2, 5, _seed = 10))
        # modules.append(TPZRandomObstruction(length, lc, circle, obstruction/2, 5))
        modules.append(TPZMultipleObstruction(length, lc, circle, obstruction/3, .3))
        modules.append(TPZNoObstruction(_length = length, _lc = lc, _module_typology = circle))
        # modules.append(TPZNoObstruction(_length = length, _lc = lc, _module_typology = circle))

    else:
        modules.append(TPZSimpleObstruction(length, lc, rec, obstruction))
        modules.append(TPZCrossObstruction(_length = length ,_lc = lc, _module_typology = rec, _radius = r, _obstruction_width = width, _obstruction_height = height))
        modules.append(TPZMultipleObstruction(length, lc, rec, obstruction/3, .3))
        modules.append(TPZSemiArcObstruction(length, lc, rec, obstruction))
        modules.append(TPZRandomObstruction(length, lc, rec, obstruction/2, 5))
        modules.append(TPZNoObstruction(_length = length, _lc = lc, _module_typology = rec))
    
    "Moving them to the right place"
    for i, module in enumerate(modules):
        module.Move(0,0, length*i)
    
    # "Joining them all in one"
    gmsh.model.occ.removeAllDuplicates()

    physical_group = [
        [(3, [1, 2]), 1, "Domain"],
        [(2, [1]), 2, "PressIn"],
        [(2, [12]), 3, "PressOut"],
        [(2, [1, 2, 3, 4, 5, 8, 9, 10, 11, 12]), 4, "NoTanVel"],
        [(2, [2, 3, 4, 5, 8, 9, 10, 11]), 5, "NoPenetration"],
        [(2, [6]), 100, "Obstruction"]
    ]

    # "Creating the elements"
    gmsh.model.occ.synchronize()

    TPZMeshModeling.CreatePhysicalGroup(physical_group)

    gmsh.model.mesh.generate(mesh_dim)
    
    TPZMeshModeling.ShowModel()

    TPZMeshModeling.WriteMeshFiles(file_name)
    TPZMeshModeling.PrintJson(json_data, file_name)

    gmsh.finalize()

if __name__ == '__main__':
    main()