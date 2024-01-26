"""
Class to create modules of length 'size'

Created by Carlos Puga: 01/13/2024
"""
#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
import gmsh

from TPZMeshModeling import TPZMeshModeling
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZModuleTypology:
    """
    Base class used to generate the module on which 
    the obstruction will be inserted. Every module is 
    created in the origin of cartesian plan

    Fields:
        - length: module length

        - lc: mesh size (gmsh requirement)

        - points: module's points (class provides it)

        - lines: module's lines (class provides it)

        - curves: module's curve loops (class provides it)

        - surfaces: module's surfaces, expcept by the one on which
         the obstruction will be inserted (class provides it)

        - obstruciton_face: module's surface on which the obstruction 
        will be inserted
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    _length: float
    _lc: float
    _points: tuple[int] = field(init=False, default_factory=tuple)
    _lines: tuple[int] = field(init=False, default_factory=tuple)
    _curves: tuple[int] = field(init=False, default_factory=tuple)
    _surfaces: list[int] = field(init=False, default_factory=list)
    _obstruction_face: int = field(init=False)

#   ****************** 
#   GETTERS & SETTERS
#   ****************** 
    @property
    def length(self)->float: return self._length
    @length.setter
    def length(self, size)->None: self._length = size

    @property
    def lc(self)->float: return self._lc
    @lc.setter
    def lc(self, LC)->None: self._lc = LC  

    @property
    def points(self)->tuple[int]: return self._points
    @points.setter
    def points(self, Points)->None: self._points = Points

    @property
    def lines(self)->tuple[int]: return self._lines
    @lines.setter
    def lines(self, Lines)->None: self._lines = Lines

    @property
    def curves(self)->tuple[int]: return self._curves
    @curves.setter
    def curves(self, Curves)->None: self._curves = Curves

    @property
    def surfaces(self)->tuple[int]: return self._surfaces
    @surfaces.setter
    def surfaces(self, Surfaces)->None: self._surfaces = Surfaces

    @property
    def obstruction_face(self)->int: return self._obstruction_face
    @obstruction_face.setter
    def obstruction_face(self, face)->None: self._obstruction_face = face    

#   ****************** 
#        METHODS
#   ******************  
    def DebugStop(self, message=''):
        raise ValueError(message + ' YOUR CHANCE TO PUT A BREAK POINT HERE')

#   ****************** 
#          BOX
#   ******************  
    def BoxPoints(self, dx: float, dy: float)->tuple[int]:
        """
        Returns the basic points of a rectangular box of dimensions '(dx, dy, length)'
        with mesh size 'lc'. 
        """
        points_coord = [ 
            [0,0,0],
            [dx, 0, 0],
            [dx, dy, 0],
            [0, dy, 0],

            [0, 0, self.length],
            [dx, 0, self.length],
            [dx, dy, self.length],
            [0, dy, self.length]
        ]

        points = TPZMeshModeling.CreatePoints(points_coord, self.lc)

        return points
    
    def BoxLines(self)->tuple[int]:
        """
        Returns the basic lines of a rectangular box using the points p1, ..., p8.
        """
        p1, p2, p3, p4, p5, p6, p7, p8 = self.points

        lines_points = [
            [p1, p2],
            [p2, p3],
            [p3, p4],
            [p4, p1],

            [p5, p6],
            [p6, p7],
            [p7, p8],
            [p8, p5],

            [p4, p8],
            [p3, p7],
            
            [p1, p5],
            [p2, p6]
        ]

        lines = TPZMeshModeling.CreateLines(lines_points)

        return lines
    
    def BoxCurveLoops(self)->tuple[int]:
        """
        Returns the basic curve loops of a rectangular box 
        using the lines l1, ..., l12
        """
        l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12 = self.lines

        curve_lines = [
            [l1, l2, l3, l4],
            [l5, l6, l7, l8],
            [l3, l10, l7, l9],
            [l1, l11, l5, l12],
            [l4, l9, l8, l11],
            [l2, l10, l6, l12]
        ] 

        curves = TPZMeshModeling.CreateCurveLoops(curve_lines)

        return curves

    def BoxSurfaces(self)->tuple[int]:
        """
        Returns the basic surfaces of rectangular box 
        usning the curve loops c1, ..., c2
        """
        c1, c2, c3, c4, c5, c6 = self.curves

        surfaces_curves = [
            [c1],
            [c2],
            [c3],
            [c4],
            [c5],
            [c6]
        ]

        surfaces = TPZMeshModeling.CreatePlanes(surfaces_curves)

        return surfaces

    def CreateBox(self, dx:float , dy: float)->None:
        """
        Create a box with dimensions 'dx', 'dy' and 'length'.
        """
        # creating the base box
        self.points = self.BoxPoints(dx, dy)
        self.lines = self.BoxLines()
        self.curves = self.BoxCurveLoops()
        s1, s2, s3, s4, s5, s6 = self.BoxSurfaces()

        self.surfaces = [s1,s3,s4,s5,s6]
        self.obstruction_face = s2

#   ****************** 
#       CYLINDER
#   ******************  
    def CylinderPoints(self, radius:float)->tuple[int]:   
        """
        Returns the cylinder points
        """
        points_coord = [
            [0, 0, 0],
            [radius, 0, 0],
            [0, radius, 0],
            [-radius, 0, 0],
            [0, -radius, 0],

            [0, 0, self.length],
            [radius, 0, self.length],
            [0, radius, self.length],
            [-radius, 0, self.length],
            [0, -radius, self.length]
        ]

        points = TPZMeshModeling.CreatePoints(points_coord, self.lc)

        return points

    def CylinderArcs(self)->tuple[int]:
        """
        Returns the lines of the cylinder
        """
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = self.points

        line_points = [
            [p2, p1, p3],
            [p3, p1, p4],
            [p4, p1, p5],
            [p5, p1, p2],

            [p7, p6, p8],
            [p8, p6, p9],
            [p9, p6, p10],
            [p10, p6, p7]
        ]

        arcs = TPZMeshModeling.CreateCircleArcs(line_points)

        gmsh.model.occ.remove([(0, p1)]) # center of the back circle
        gmsh.model.occ.remove([(0, p6)]) # center of the front circle
        # this was made due to some problems we faced creating the mesh

        return arcs
    
    def CylinderCurves(self)->tuple[int]:
        """
        Returns the cylinder curve loops
        """
        l1, l2, l3, l4, l5, l6, l7, l8 = self.lines

        back = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
        front = gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])

        curves = (back, front)

        return curves
    
    def CylinderSurfaces(self)->list[int]:
        """
        Returns the cylinder surfaces
        """
        c1, c2 = self.curves

        back = gmsh.model.occ.addPlaneSurface([c1])
        front = gmsh.model.occ.addPlaneSurface([c2])

        contour = gmsh.model.occ.addThruSections([back, front], makeSolid=False)

        surfaces = [back, front] + [c[1] for c in contour]

        return surfaces


    def CreateCylinder(self, radius:float)->None:
        """
        Create a cylinder with 'radius'
        """

        self.points = self.CylinderPoints(radius)
        self.lines = self.CylinderArcs()
        self.curves = self.CylinderCurves()
        s1, s2, s3, s4, s5, s6 = self.CylinderSurfaces()

        self.surfaces = [s1, s3, s4, s5, s6]
        self.obstruction_face = s2
# %%
