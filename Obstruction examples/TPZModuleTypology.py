"""
Class to create modules of length 'size'

Created by Carlos Puga: 01/13/2024
"""
#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
import gmsh
import sys
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
        back_lower_left = gmsh.model.occ.addPoint(0, 0, 0, self.lc)
        back_lower_right = gmsh.model.occ.addPoint(dx, 0, 0, self.lc)
        back_upper_right = gmsh.model.occ.addPoint(dx, dy, 0, self.lc)
        back_upper_left = gmsh.model.occ.addPoint(0, dy, 0, self.lc)

        front_lower_left = gmsh.model.occ.addPoint(0, 0, self.length, self.lc)
        front_lower_right = gmsh.model.occ.addPoint(dx, 0, self.length, self.lc)
        front_upper_right = gmsh.model.occ.addPoint(dx, dy, self.length, self.lc)
        front_upper_left = gmsh.model.occ.addPoint(0, dy, self.length, self.lc)

        points = (back_lower_left, back_lower_right, back_upper_right, back_upper_left, front_lower_left, front_lower_right, front_upper_right, front_upper_left)

        return points
    
    def BoxLines(self)->tuple[int]:
        """
        Returns the basic lines of a rectangular box using the points p1, ..., p8.
        """
        p1, p2, p3, p4, p5, p6, p7, p8 = self.points

        back_lower = gmsh.model.occ.addLine(p1, p2)
        back_right = gmsh.model.occ.addLine(p2, p3)
        back_upper = gmsh.model.occ.addLine(p3, p4)
        back_left = gmsh.model.occ.addLine(p4, p1)

        front_lower = gmsh.model.occ.addLine(p5, p6)
        front_right = gmsh.model.occ.addLine(p6, p7)
        front_upper = gmsh.model.occ.addLine(p7, p8)
        front_left = gmsh.model.occ.addLine(p8, p5)

        upper_left = gmsh.model.occ.addLine(p4, p8)
        upper_right = gmsh.model.occ.addLine(p3, p7)

        lower_left = gmsh.model.occ.addLine(p1, p5)
        lower_right = gmsh.model.occ.addLine(p2, p6)

        lines = (back_lower, back_right, back_upper, back_left, front_lower, front_right, front_upper, front_left, upper_left, upper_right, lower_left, lower_right)

        return lines
    
    def BoxCurveLoops(self)->tuple[int]:
        """
        Returns the basic curve loops of a rectangular box 
        using the lines l1, ..., l12
        """
        l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12 = self.lines

        back = gmsh.model.occ.addCurveLoop([l1,l2,l3,l4])
        front = gmsh.model.occ.addCurveLoop([l5,l6,l7,l8])

        upper = gmsh.model.occ.addCurveLoop([l3,l10,l7,l9])
        lower = gmsh.model.occ.addCurveLoop([l1, l11, l5, l12])

        left = gmsh.model.occ.addCurveLoop([l4, l9,l8, l11])
        right = gmsh.model.occ.addCurveLoop([l2,l10,l6,l12])

        curves = (back, front, upper, lower, left, right)

        return curves

    def BoxSurfaces(self)->tuple[int]:
        """
        Returns the basic surfaces of rectangular box 
        usning the curve loops c1, ..., c2
        """
        c1, c2, c3, c4, c5, c6 = self.curves

        back = gmsh.model.occ.addPlaneSurface([c1])
        front = gmsh.model.occ.addPlaneSurface([c2])

        lower = gmsh.model.occ.addPlaneSurface([c3])
        upper = gmsh.model.occ.addPlaneSurface([c4])

        left = gmsh.model.occ.addPlaneSurface([c5])
        right = gmsh.model.occ.addPlaneSurface([c6])

        surfaces = (back, front, lower, upper, left, right)

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
        back_center = gmsh.model.occ.addPoint(0, 0, 0, self.lc)
        back_right = gmsh.model.occ.addPoint(radius, 0, 0, self.lc)
        back_upper = gmsh.model.occ.addPoint(0, radius, 0, self.lc)
        back_left = gmsh.model.occ.addPoint(-radius, 0, 0, self.lc)
        back_lower = gmsh.model.occ.addPoint(0, -radius, 0, self.lc)

        front_center = gmsh.model.occ.addPoint(0, 0, self.length, self.lc)
        front_right = gmsh.model.occ.addPoint(radius, 0, self.length, self.lc)
        front_upper = gmsh.model.occ.addPoint(0, radius, self.length, self.lc)
        front_left = gmsh.model.occ.addPoint(-radius, 0, self.length, self.lc)
        front_lower = gmsh.model.occ.addPoint(0, -radius, self.length, self.lc)

        points = (back_center, back_right, back_upper, back_left, back_lower, front_center, front_right, front_upper, front_left, front_lower)

        return points

    def CylinderArcs(self)->tuple[int]:
        """
        Returns the lines of the cylinder
        """
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = self.points

        back_upper_right = gmsh.model.occ.addCircleArc(p2, p1, p3)
        back_upper_left = gmsh.model.occ.addCircleArc(p3, p1, p4)
        back_lower_left = gmsh.model.occ.addCircleArc(p4, p1, p5)
        back_lower_right = gmsh.model.occ.addCircleArc(p5, p1, p2)

        front_upper_right = gmsh.model.occ.addCircleArc(p7, p6, p8)
        front_upper_left = gmsh.model.occ.addCircleArc(p8, p6, p9)
        front_lower_left = gmsh.model.occ.addCircleArc(p9, p6, p10)
        front_lower_right = gmsh.model.occ.addCircleArc(p10, p6, p7)

        gmsh.model.occ.remove([(0, p1)])
        gmsh.model.occ.remove([(0, p6)])

        lines = (back_upper_right, back_upper_left, back_lower_left, back_lower_right, front_upper_right, front_upper_left, front_lower_left, front_lower_right)

        return lines
    
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
