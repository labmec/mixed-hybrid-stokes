"""
Code to generate a mesh of a domain with fractures
"""
# -----------------------
#   Importing modules
# -----------------------
from TPZMeshModeling import TPZMeshModeling
import gmsh

# -----------------------
#       Functions
# -----------------------
def DispX0(x0:tuple[float], disp:tuple[float])->list[float]:
    """
    Return the point x0 displaced by the vector disp
    """
    return [X + d for X, d in zip(x0, disp)]

# -----------------------
#   Main function
# -----------------------
def main()->None:
    # initial coordinate
    x0 = [0., 0., 0.]

    # length values
    l1:float = 10
    l2:float = l1 / 2 
    l3:float = l2 / 3
    l4:float = l3 * 1.05

    # point coordinates
    x1: float = l1 - l4 - 1 
    x2: float = l1 / 2 - 0.5
    x3: float = x2 - 0.5
    x4: float = l1 - l4 + 0.3

    y1: float = 2 / 3 * l2 + 0.3
    y2: float = 1 / 3 * l2 + 0.1
    y3: float = 1 / 3 * l2 - 0.1
    y4: float = 1 / 3  * l2 * 0.8

    z = 0

    # mesh size
    lc = 5e-1

    # show mesh bool 
    # change this to hide or show the mesh 
    show = True

    TPZMeshModeling.Begin()
    TPZMeshModeling.TurnOnLabels('point', 'line', 'surface')
    TPZMeshModeling.TurnOnRepresentation('surfaces')

    # Creating the surface domain
    # point coordinates
    M1_coord = [
        x0,
        DispX0(x0, (l1, 0., 0.)),
        DispX0(x0, (l1, (l2 - l3)/2, 0.)),
        DispX0(x0, (l1, (l2 + l3)/2, 0.)),
        DispX0(x0, (l1, l2, 0.)),
        DispX0(x0, (0, l2, 0.)),
    ]

    M2_coord = [
        DispX0(M1_coord[2], (-l4, 0., 0.),),
        DispX0(M1_coord[2], (-l4, l3, 0.),)
    ]

    point_cord = M1_coord + M2_coord

    p1, p2, p3, p4, p5, p6, p7, p8 = TPZMeshModeling.CreatePoints(point_cord, lc)

    # lines
    line_points = [
        [p1, p2],
        [p2, p3],
        [p3, p7],
        [p7, p8],
        [p8, p4],
        [p4, p5],
        [p5, p6],
        [p6, p1],

        [p3, p4],
    ]

    l1, l2, l3, l4, l5, l6, l7, l8, l9 = TPZMeshModeling.CreateLines(line_points)

    # curve loops
    c1 = TPZMeshModeling.CreateCurveLoops([[l1, l2, l3, l4, l5, l6, l7, l8]])
    c2 = TPZMeshModeling.CreateCurveLoops([[l3, l4, l5, l9]])

    # plane 
    plane1 = TPZMeshModeling.CreatePlanes([c1])[0]
    plane2 = TPZMeshModeling.CreatePlanes([c2])[0]

    # fracture points
    frac_points = [
        (x1, y1, z),
        (x2, y2, z),
        (x3, y3, z),
        (x4, y4, z)
    ]

    fp1, fp2, fp3, fp4 = TPZMeshModeling.CreatePoints(frac_points, lc)

    # fracture lines
    fl1 = TPZMeshModeling.CreateLines([[fp1, fp2]])[0]
    fl2 = TPZMeshModeling.CreateLines([[fp3, fp4]])[0]

    TPZMeshModeling.Synchronize()

    gmsh.model.mesh.embed(1, [fl1, fl2], 2, plane1)

    TPZMeshModeling.Synchronize()

    pg = [
        [(2, [plane1]), 1, "Domain1"],
        [(2, [plane2]), 2, "Domain2"],
        [(1, [l2]), 3, "FlowOut"],
        [(1, [l6]), 4, "FlowIn"],
        [(1, [l1, l9, l7, l8]), 5, "NoFlow"],
        [(1, [fl1]), 6, "Fracture1"],
        [(1, [fl2]), 7, "Fracture2"]
    ]

    TPZMeshModeling.CreatePhysicalGroup(pg)

    TPZMeshModeling.CreateMesh(2)

    if show:
        TPZMeshModeling.ShowModel()

    TPZMeshModeling.WriteMeshFiles("FractureExample", ".msh")

    TPZMeshModeling.End()


if __name__ == "__main__":
    main()