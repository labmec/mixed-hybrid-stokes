"""
Class to create obstruction with semi arcs on 
the module's front face

Created by Carlos Puga: 01/15/2024
"""

#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
from typing import ClassVar
import numpy as np
import gmsh
import sys

from TPZModuleTypology import TPZModuleTypology
from TPZMeshModeling import TPZMeshModeling
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZSemiArcObstruction(TPZModuleTypology):
    """
    This class creates a obstruction with semi arcs in the 
    middle of the front (positive direction of z axis) face.
    
    It inherits a base creator of typologies (TPZModuleTypology)
    to generates the module body and insert the obstruction.

    Fields: 
        - module_typology: information about the type 
        of module to be constructed. Currently is only availabe rectangular 
        modules. More information are required: 
            * Rectangular -> {dx: length on x-direction, dy: length on y-direction}
            * Circular -> {radius: modules'radius}

        - obstruction_radius: distance from the center of the module to the semi arcs

        - obstruction_cx: x coordinate to place the obstruction

        - obstruction_cy: y coordinate to place the obstruction

        - id: obstruction volume's identification (provided by the class itself)
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    multiplicative_constant:ClassVar = 1.2

    _module_typology: tuple[str, dict]
    _obstruction_radius: float
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

            if 2*self.obstruction_radius > dx or 2*self.obstruction_radius > dy:
                self.DebugStop('ERROR: obstruction radius not compatible with box dimensions!')

            self.CreateBox(dx, dy)

        elif m_type == 'Circular':
            if not 'radius' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            radius = m_characteristics['radius']

            self.obstruction_cx = 0
            self.obstruction_cy = 0

            if self.obstruction_radius > radius:
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
        new_surfaces = gmsh.model.occ.fragment([(2, self.obstruction_face)], [(2, obs) for obs in obstruction])

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
        r = self.obstruction_radius
        l = self.length
        lc = self.lc
        cons = self.multiplicative_constant

        dx = r*np.cos(np.deg2rad(60))
        dy = r*np.sin(np.deg2rad(60))

        points_coord = [
            [cx, cy, l],

            [cx + dx, cy + dy, l],
            [cx + cons * dx, cy + cons * dy, l],
            [cx - cons * dx, cy + cons * dy, l],
            [cx - dx, cy + dy, l],
            
            [cx + dx, cy - dy, l],
            [cx + cons * dx, cy - cons * dy, l],
            [cx - dx, cy - dy, l],
            [cx - cons * dx, cy - cons * dy, l]
        ]

        dx = r*np.cos(np.deg2rad(30))
        dy = r*np.sin(np.deg2rad(30))

        points_coord += [
            [cx + dx, cy + dy, l],
            [cx + cons * dx, cy + cons * dy, l],
            [cx + dx, cy - dy, l],
            [cx + cons * dx, cy - cons * dy, l],

            [cx - dx, cy + dy, l],
            [cx - cons * dx, cy + cons * dy, l],
            [cx - dx, cy - dy, l],
            [cx - cons * dx, cy - cons * dy, l]
        ]


        ob_points = TPZMeshModeling.CreatePoints(points_coord, lc)

        return ob_points

    def ObstructionArcs(self, points: list[int])->list[int]:
        """
        Returns the lines belonging to the obstruction
        """
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17 = points

        line_points = [
            [[p2, p3], [p3, p1, p4], [p4, p5], [p5, p1, p2]],
            [[p7, p6], [p6, p1, p8], [p8, p9], [p9, p1, p7]],
            [[p12, p13], [p13, p1, p11], [p11, p10], [p10, p1, p12]],
            [[p14, p15], [p15, p1, p17], [p17, p16], [p16, p1, p14]]
        ]

        ob_arcs = []
        for group in line_points:
            g1, g2, g3, g4 = group
            
            arc = []
            arc.append(gmsh.model.occ.addLine(g1[0], g1[1]))
            arc.append(gmsh.model.occ.addCircleArc(g2[0], g2[1], g2[2]))
            arc.append(gmsh.model.occ.addLine(g3[0], g3[1]))
            arc.append(gmsh.model.occ.addCircleArc(g4[0], g4[1], g4[2]))

            ob_arcs.append(arc)

        gmsh.model.occ.remove([(0, p1)])

        return ob_arcs

    def CreateObstruction(self)->int:
        """
        Returns the obstruction surface id
        """
        ob_points = self.ObstructionPoints()
        ob_arcs = self.ObstructionArcs(ob_points)
        
        obstruction = []
        for arc in ob_arcs:
            curve = gmsh.model.occ.addCurveLoop(arc)
        
            obstruction_surface = gmsh.model.occ.addPlaneSurface([curve])
            obstruction.append(obstruction_surface)

        return obstruction
    
    def Move(self, dx:float, dy:float, dz:float):
        """
        Moves the module (dx, dy, dz)
        """
        gmsh.model.occ.translate([ (3,self.id)], dx, dy, dz)
