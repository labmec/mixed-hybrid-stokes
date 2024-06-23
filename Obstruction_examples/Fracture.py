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
    # ATTENTION: remeber to always create the mesh during the show model 
    # since now the model is shown three times, it is possible that the mesh is not created in one of them. 
    # Just type 2 on your keyboard to be sure that the mesh is created

    # initial coordinate
    x0 = [0., 0., 0.]

    # length values
    l1:float = 10 # domain 1 length 
    l2:float = l1 / 2 # domain 1 height
    l3:float = l2 / 3 # domain 2 height
    l4:float = l3 * 1.05 # domain 2 length

    # point coordinates
    x1: float = l1 - l4 - 1  # fracture 1 first x coordinate
    x2: float = l1 / 2 - 0.5 # fracture 1 second x coordinate
    x3: float = x2 - 0.5     # fracture 2 first x coordinate
    x4: float = l1 - l4 + 0.3 # fracture 2 second x coordinate

    y1: float = 2 / 3 * l2 + 0.3 # fracture 1 first y coordinate
    y2: float = 1 / 3 * l2 + 0.1 # fracture 1 second y coordinate
    y3: float = 1 / 3 * l2 - 0.1 # fracture 2 first y coordinate
    y4: float = 1 / 3  * l2 * 0.8 # fracture 2 second y coordinate

    z = 0 # z coordinate

    # mesh size
    lc = 5e-1

    # show mesh bool 
    # change this to hide or show the mesh 
    show = True

    TPZMeshModeling.Begin()
    TPZMeshModeling.TurnOnLabels('point', 'line', 'surface') # turn on the labels of the points, lines and surfaces
    TPZMeshModeling.TurnOnRepresentation('surfaces') # turn on the geometrical representation of the surfaces

    # Creating the surface domain
    # point coordinates
    M1_coord = [
        x0, # p1
        DispX0(x0, (l1, 0., 0.)), # p2
        DispX0(x0, (l1, (l2 - l3)/2, 0.)), # p3
        DispX0(x0, (l1, (l2 + l3)/2, 0.)), # p4
        DispX0(x0, (l1, l2, 0.)), # p5
        DispX0(x0, (0, l2, 0.)), # p6
    ]

    M2_coord = [
        DispX0(M1_coord[2], (-l4, 0., 0.),), # p7
        DispX0(M1_coord[2], (-l4, l3, 0.),) # p8
    ]

    point_cord = M1_coord + M2_coord

    p1, p2, p3, p4, p5, p6, p7, p8 = TPZMeshModeling.CreatePoints(point_cord, lc)

    # lines
    line_points = [
        [p1, p2], # l1 
        [p2, p3], # l2
        [p3, p7], # l3
        [p7, p8], # l4
        [p8, p4], # l5 
        [p4, p5], # l6
        [p5, p6], # l7
        [p6, p1], # l8

        [p3, p4], # l9
    ]

    l1, l2, l3, l4, l5, l6, l7, l8, l9 = TPZMeshModeling.CreateLines(line_points)

    # curve loops
    c1 = TPZMeshModeling.CreateCurveLoops([[l1, l2, l3, l4, l5, l6, l7, l8]]) # domain 1 curve loop 
    c2 = TPZMeshModeling.CreateCurveLoops([[l3, l4, l5, l9]]) # domain 2 curve loop

    # plane 
    plane1 = TPZMeshModeling.CreatePlanes([c1])[0] # domain 1 plane surface
    plane2 = TPZMeshModeling.CreatePlanes([c2])[0] # domain 2 plane surface

    # fracture points
    frac_points = [
        (x1, y1, z), # p9
        (x2, y2, z), # p10
        (x3, y3, z), # p11
        (x4, y4, z) # p12
    ]

    fp1, fp2, fp3, fp4 = TPZMeshModeling.CreatePoints(frac_points, lc)

    # fracture lines
    fl1 = TPZMeshModeling.CreateLines([[fp1, fp2]])[0] # fracture 1 line
    fl2 = TPZMeshModeling.CreateLines([[fp3, fp4]])[0] # fracture 2 line

    TPZMeshModeling.Synchronize()

    gmsh.model.mesh.embed(1, [fl1, fl2], 2, plane1) # embed the fracture lines in the domain 1 plane

    TPZMeshModeling.Synchronize()

    if show:
        TPZMeshModeling.ShowModel() # showing the model without the physical groups

    # physical groups
    pg = [
        [(2, [plane1]), 1, "Domain1"],
        [(2, [plane2]), 2, "Domain2"],
        [(1, [l2]), 3, "FlowOut"],
        [(1, [l6]), 4, "FlowIn"],
        [(1, [l1, l9, l7, l8]), 5, "NoFlow"],
        [(1, [fl1]), 6, "Fracture1"],
        [(1, [fl2]), 7, "Fracture2"]
    ]

    TPZMeshModeling.CreatePhysicalGroup(pg) # creating the physical groups

    TPZMeshModeling.CreateMesh(2) # creating the mesh

    if show:
        # showing the model only with the physical group tags 
        TPZMeshModeling.TurnOffLabels('point')
        TPZMeshModeling.SetDescription("physical_tag")
        TPZMeshModeling.ShowModel()

        # showing the model only with the physical group names
        TPZMeshModeling.SetDescription("physical_name")
        TPZMeshModeling.ShowModel()

    TPZMeshModeling.WriteMeshFiles("FractureExample", ".msh")

    TPZMeshModeling.End()


if __name__ == "__main__":
    main()