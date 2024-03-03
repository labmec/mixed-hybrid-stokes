"""
Class to create several circular obstructions on 
the module's front face

Created by Carlos Puga: 01/15/2024
"""

#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
import numpy as np
import gmsh

from TPZModuleTypology import TPZModuleTypology
from TPZMeshModeling import TPZMeshModeling
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZMultipleObstruction(TPZModuleTypology):
    """
    This class creates several circular obstructions in the 
    middle of the front (positive direction of z axis) face.
    
    It inherits a base creator of typologies (TPZModuleTypology)
    to generates the module body and insert the obstruction.

    Fields: 
        - module_typology: information about the type 
        of module to be constructed. Currently is only availabe rectangular 
        modules. More information are required: 
            * Rectangular -> {dx: length on x-direction, dy: length on y-direction}
            * Circular -> {radius: modules'radius}

        - obstruction_radius: size of the obstrcution radius

        - obstruction_distance: distance from the center of the domain

        - obstruction_cx: x coordinate to place the obstruction

        - obstruction_cy: y coordinate to place the obstruction

        - id: obstruction volume's identification (provided by the class itself)
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    _module_typology: tuple[str, dict]
    _obstruction_radius: float
    _obstruction_distance: float
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
    def obstruction_radius(self)->float: return self._obstruction_radius
    @obstruction_radius.setter
    def obstruction_radius(self, R)->None: self._obstruction_radius = R

    @property
    def obstruction_distance(self)->float: return self._obstruction_distance
    @obstruction_distance.setter
    def obstruction_distance(self, d)->None: self._obstruction_distance = d

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
        m_type, m_characteristics = self.module_typology
        keys = list(m_characteristics.keys())

        if m_type == 'Rectangular':
            if not 'dx' in keys and not 'dy' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            dx = m_characteristics['dx']
            dy = m_characteristics['dy']

            self.obstruction_cx = dx/2
            self.obstruction_cy = dy/2

            if (self.obstruction_radius+self.obstruction_distance) > dx or (self.obstruction_radius+self.obstruction_distance) > dy:
                self.DebugStop('ERROR: obstruction radius not compatible with box dimensions!')

            self.CreateBox(dx, dy)

        elif m_type == 'Circular':
            if not 'radius' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            radius = m_characteristics['radius']

            self.obstruction_cx = 0
            self.obstruction_cy = 0

            if (self.obstruction_radius+self.obstruction_distance) > radius:
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
        obstructions = self.CreateObstruction()

        # calculating the boolean difference between the domain and the obstruction faces  
        new_surfaces = gmsh.model.occ.fragment([(2, self.obstruction_face)], [(2, obs) for obs in obstructions])

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
        cx = self.obstruction_cx
        cy = self.obstruction_cy
        r = self.obstruction_radius
        l = self.length
        lc = self.lc

        point_coords = [
            [cx, cy, l],
            [cx + r, cy, l],
            [cx, cy + r, l],
            [cx - r, cy, l],
            [cx, cy - r, l]
        ]

        ob_points = TPZMeshModeling.CreatePoints(point_coords, lc)

        return ob_points

    def ObstructionArcs(self, points: list[int])->tuple[int]:
        """
        Returns the lines belonging to the obstruction
        """
        p1, p2, p3, p4, p5 = points

        arc_points = [
            [p2, p1, p3],
            [p3, p1, p4],
            [p4, p1, p5],
            [p5, p1, p2]
        ]

        ob_arcs = TPZMeshModeling.CreateCircleArcs(arc_points)

        gmsh.model.occ.remove([(0, p1)])

        return ob_arcs

    def CreateObstruction(self)->list[int]:
        """
        Returns the obstruction surface id
        """
        angles = [0,45,90,135,180,225,270,315]
        
        origin_x = self.obstruction_cx
        origin_y = self.obstruction_cy
        r = self.obstruction_distance
        
        multiple_obs = []

        center_points = self.ObstructionPoints()
        center_arcs = self.ObstructionArcs(center_points)
        center_curves = gmsh.model.occ.addCurveLoop(center_arcs)        
        center_surface = gmsh.model.occ.addPlaneSurface([center_curves])

        multiple_obs.append(center_surface)

        for theta in angles:
            self.obstruction_cx = r*np.cos(np.deg2rad(theta)) + origin_x
            self.obstruction_cy = r*np.sin(np.deg2rad(theta)) + origin_y

            points = self.ObstructionPoints()
            arcs = self.ObstructionArcs(points) 
            curves = gmsh.model.occ.addCurveLoop(arcs)
            circle_surface = gmsh.model.occ.addPlaneSurface([curves])

            multiple_obs.append(circle_surface)

        return multiple_obs
    
    def Move(self, dx:float, dy:float, dz:float):
        """
        Moves the module (dx, dy, dz)
        """
        gmsh.model.occ.translate([ (3,self.id)], dx, dy, dz)