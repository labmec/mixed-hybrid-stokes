"""
Class to create a cross obstruction on 
the module's front face

Created by Carlos Puga: 01/15/2024
"""

#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
import gmsh

from TPZModuleTypology import TPZModuleTypology
from TPZMeshModeling import TPZMeshModeling
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZCrossObstruction(TPZModuleTypology):
    """
    This class creates a cross obstruction in the 
    middle of the front (positive direction of z axis) face.
    
    It inherits a base creator of typologies (TPZModuleTypology)
    to generates the module body and insert the obstruction.

    Fields: 
        - module_typology: information about the type 
        of module to be constructed. Currently is only availabe rectangular 
        modules. More information are required: 
            * Rectangular -> {dx: length on x-direction, dy: length on y-direction}

        - obstruction_width: obstruction width

        - obstruction_height: obstruction height

        - radius: cross curved edge's radius 

        - obstruction_cx: x coordinate to place the obstruction

        - obstruction_cy: y coordinate to place the obstruction

        - id: obstruction volume's identification (provided by the class itself)
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    _module_typology: tuple[str, dict]
    _obstruction_width: float
    _obstruction_height: float
    _radius: float
    _obstruction_cx: float = field(init=False)
    _obstruction_cy: float = field(init=False)
    _id: int = field(init=False)

    def __post_init__(self):
        self.CreateDomain()
    
#   ****************** 
#   GETTERS & SETTERS
#   ******************  
    @property
    def module_typology(self)->tuple[str, dict]: return self._module_typology
    @module_typology.setter
    def module_typology(self, Module)->None: self._module_typology = Module

    @property
    def obstruction_width(self)->float: return self._obstruction_width
    @obstruction_width.setter
    def obstruction_width(self, W)->None: self._obstruction_width = W

    @property
    def obstruction_height(self)->float: return self._obstruction_height
    @obstruction_height.setter
    def obstruction_height(self, H)->None: self._obstruction_height = H

    @property
    def radius(self)->float: return self._radius
    @radius.setter
    def radius(self, R)->None: self._radius = R

    @property
    def obstruction_cx(self)->float: return self._obstruction_cx
    @obstruction_cx.setter
    def obstruction_cx(self, cx)->None: self._obstruction_cx = cx

    @property
    def obstruction_cy(self)->float: return self._obstruction_cy
    @obstruction_cy.setter
    def obstruction_cy(self, cy)->None: self._obstruction_cy = cy

    @property
    def id(self)->int: return self._id
    @id.setter
    def id(self, ID)->None: self._id = ID

#   ****************** 
#        METHODS
#   ******************  
    def CheckTypology(self)->None:
        """
        Checks whether the module is rectangular or circular 
        """
        m_type, m_characteristics = self.module_typology
        keys = list(m_characteristics.keys())

        if m_type == 'Rectangular':
            if not 'dx' in keys and not 'dy' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            dx = m_characteristics['dx']
            dy = m_characteristics['dy']

            self.obstruction_cx = dx/2
            self.obstruction_cy = dy/2

            if (self.obstruction_width + self.radius) > dx/2 or (self.obstruction_height + self.radius) > dy/2:
                self.DebugStop('ERROR: obstruction radius not compatible with box dimensions!')

            self.CreateBox(dx, dy)

        elif m_type == 'Circular':
            if not 'radius' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            radius = m_characteristics['radius']

            self.obstruction_cx = 0
            self.obstruction_cy = 0

            if (self.obstruction_width + self.radius) > radius or (self.obstruction_height + self.radius) > radius:
                self.DebugStop('ERROR: obstruction radius not compatible with cylinder dimensions!')

            self.CreateCylinder(radius)

        else:
            self.DebugStop(f"ERROR: I'm sorry the typology {m_type} has not been implemented yet")


    def CreateDomain(self)->None:
        """
        Generates the module with a circular obstruction on gmsh.
        """
        self.CheckTypology()

        # creating the obstruction on surface obstruction_face
        obstruction = self.CreateObstruction()

        # calculating the boolean difference between the domain and the obstruction faces  
        new_surfaces = gmsh.model.occ.fragment([(2, self.obstruction_face)], [(2, obstruction)])

        domain_surfaces = self.surfaces + [surface[1] for surface in new_surfaces[0]]

        # creating the module volume
        surface_loop = gmsh.model.occ.addSurfaceLoop(domain_surfaces)
        domain_volume = gmsh.model.occ.addVolume([surface_loop])

        self.id = domain_volume

#   ****************** 
#        OBSTRUCTION
#   ******************  
    def ObstructionPoints(self)->tuple[int]:
        """
        Returns the points belonging to the obstruction
        """
        cx, cy = self.obstruction_cx, self.obstruction_cy
        dx, dy = self.obstruction_width, self.obstruction_height
        l = self.length
        r = self.radius
        lc = self.lc

        point_coord = [
            [cx + dx, cy, l],
            [cx + dx, cy - r, l],
            [cx + dx + r, cy, l],
            [cx + dx, cy + r, l],

            [cx, cy + dy, l],
            [cx + r, cy + dy, l],
            [cx, cy + dy + r, l],
            [cx - r, cy + dy, l],

            [cx - dx, cy, l],
            [cx - dx, cy + r, l],
            [cx - dx - r, cy, l],
            [cx - dx, cy - r, l],

            [cx, cy - dy, l],
            [cx - r, cy - dy, l],
            [cx, cy - dy - r, l],
            [cx + r, cy - dy, l],

            [cx + r, cy + r, l],
            [cx - r, cy + r, l],
            [cx - r, cy - r, l],
            [cx + r, cy - r, l]
        ]

        ob_points = TPZMeshModeling.CreatePoints(point_coord, lc)

        return ob_points

    def ObstructionArcs(self, points:list[int])->list[int]:
        """
        Returns the lines belonging to the obstruction
        """
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20 = points

        points = [
            {"arcs": [[p2, p1, p3], [p3, p1, p4]], 
             "lines": [[p4, p17], [p17, p6]]},

            {"arcs":[[p6, p5, p7], [p7, p5, p8]],
             "lines": [[p8, p18], [p18, p10]]},

            {"arcs": [[p10, p9, p11], [p11, p9, p12]],
             "lines": [[p12, p19], [p19, p14]]},

            {"arcs": [[p14, p13, p15], [p15, p13, p16]],
             "lines": [[p16, p20], [p20, p2]]} 
        ]

        obs_lines = []
        for group in points:
            gp_arcs = group["arcs"]
            gp_lines = group["lines"]

            a = TPZMeshModeling.CreateCircleArcs(gp_arcs)
            l = TPZMeshModeling.CreateLines(gp_lines)
            
            obs_lines += a
            obs_lines += l

        gmsh.model.occ.remove([(0, p1)])
        gmsh.model.occ.remove([(0, p5)])
        gmsh.model.occ.remove([(0, p9)])
        gmsh.model.occ.remove([(0, p13)])

        return obs_lines

    def CreateObstruction(self)->int:
        """
        Returns the obstruction surface id
        """
        points = self.ObstructionPoints()
        curves = self.ObstructionArcs(points)        
        c1 = gmsh.model.occ.addCurveLoop(curves)
        
        obstruction_surface = gmsh.model.occ.addPlaneSurface([c1])

        return obstruction_surface
    
    def Move(self, dx:float, dy:float, dz:float):
        """
        Moves the module (dx, dy, dz)
        """
        gmsh.model.occ.translate([ (3,self.id)], dx, dy, dz)