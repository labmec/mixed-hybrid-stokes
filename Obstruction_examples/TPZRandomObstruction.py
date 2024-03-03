"""
Class to create random circular obstructions on 
the module's front face

Created by Carlos Puga: 01/15/2024
"""

#%% ****************** 
#   IMPORTED MODULES
#   ******************
import gmsh
import random
import numpy as np
from datetime import datetime
from dataclasses import dataclass, field

from TPZModuleTypology import TPZModuleTypology
from TPZMeshModeling import TPZMeshModeling
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZRandomObstruction(TPZModuleTypology):
    """
    This class creates a random circular obstructions in the 
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

        - number_of_obstrucitons: maximum number of obstructions

        - obstruction_cx: x coordinate to place the obstruction

        - obstruction_cy: y coordinate to place the obstruction

        - id: obstruction volume's identification (provided by the class itself)
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    _module_typology: tuple[str, dict]
    _obstruction_radius: float
    _number_of_obstructions: int
    _seed: int = None
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
    def number_of_obstructions(self)->int: return self._number_of_obstructions
    @number_of_obstructions.setter 
    def number_of_obstructions(self, number)->None: self._number_of_obstructions = number

    @property
    def seed(self)->int: return self._seed
    @seed.setter
    def seed(self, s)->None: self._seed = s

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
                self.DebugStop('Not enough information to create the box!')

            dx = m_characteristics['dx']
            dy = m_characteristics['dy']

            self.obstruction_cx = dx/2
            self.obstruction_cy = dy/2

            self.CreateBox(dx, dy)

        elif m_type == 'Circular':
            if not 'radius' in keys:
                self.DebugStop('Not enough information to create the cylinder!')

            radius = m_characteristics['radius']

            self.obstruction_cx = 0
            self.obstruction_cy = 0

            self.CreateCylinder(radius)

        else:
            self.DebugStop(f"I'm sorry the typology {m_type} has not been implemented yet!")

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

    def EuclideanDistance(self, xa: float, ya: float, xb: float, yb: float)->float:
        """
        Returns the euclidean distance between the points (xa, ya) and (xb, yb)
        """
        Xa = np.array([xa,ya])
        Xb = np.array([xb,yb])

        return np.linalg.norm(Xa-Xb)

    def GetObstructionDomain(self, module_type:str, module_char:dict[str,float])->tuple[float]:
        """
        Returns the domain range in which the obstructions can be inserted
        """
        if module_type == 'Rectangular':
            dx = module_char['dx']
            dy = module_char['dy']

            domX = .9*(dx-self.obstruction_cx-self.obstruction_radius)
            domY = .9*(dy-self.obstruction_cy-self.obstruction_radius)

        elif module_type == 'Circular':
            radius = module_char['radius']

            domX =  .85*(radius-self.obstruction_radius)
            domY = .85*(radius-self.obstruction_radius)
        else:
            self.DebugStop(f"I'm sorry, the typology {module_type} has not been implemented yet!")

        return domX, domY

    def NoOverlappingCircles(self):
        """
        Returns a list containing the (x,y) coordinates of 'n_samples' random circles, with no overlapping between them. The range of the coordinates is within the ['domainX' x 'domainY'] values and the circles have radius = 'radius'. In case of cylindrical domains, it is also verified whether the coordinates fit the domain or not.

        Return: circleList
        """
        counter = 0
        circleList = []

        module_type, module_char = self.module_typology

        domX, domY = self.GetObstructionDomain(module_type, module_char)

        if self.seed:
            random.seed(self.seed)
        else:
            random.seed(datetime.now().timestamp())

        while len(circleList) < self.number_of_obstructions and counter < 1000:
            counter += 1

            x_mult = (-1)**random.randint(0,1)
            y_mult = (-1)**random.randint(0,1)

            x = (x_mult)*random.uniform(0, domX)
            y = (y_mult)*random.uniform(0, domY)

            if module_type == 'Rectangular':
                if not any((Xcenter, Ycenter) for Xcenter, Ycenter in circleList if self.EuclideanDistance(x, y, Xcenter, Ycenter)< 2.5*self.obstruction_radius):
                        circleList.append((x, y)) 
            
            elif module_type == 'Circular':
                if not any((Xcenter, Ycenter) for Xcenter, Ycenter in circleList if self.EuclideanDistance(x, y, Xcenter, Ycenter)< 2.5*self.obstruction_radius):
                    radius = module_char['radius']
                    if (x)**2 + (y)**2 < (.75*radius)**2:
                        circleList.append((x, y))  

        return circleList

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

        point_coord = [
            [cx, cy, l],
            [cx + r, cy, l],
            [cx, cy + r, l],
            [cx - r, cy, l],
            [cx, cy - r, l]
        ]

        ob_points = TPZMeshModeling.CreatePoints(point_coord, lc)

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

    def CreateObstruction(self)->int:
        """
        Returns the obstruction surface id
        """
        origin_x = self.obstruction_cx
        origin_y = self.obstruction_cy
        
        obstruction_coordinates = self.NoOverlappingCircles()

        obstructions = []

        for coordinates in obstruction_coordinates:
            dx, dy = coordinates

            self.obstruction_cx = origin_x + dx
            self.obstruction_cy = origin_y + dy

            ob_points = self.ObstructionPoints()
            ob_arc = self.ObstructionArcs(ob_points)
            
            curve = gmsh.model.occ.addCurveLoop(ob_arc)
            
            obstruction_surface = gmsh.model.occ.addPlaneSurface([curve])

            obstructions.append(obstruction_surface)

        return obstructions
    
    def Move(self, dx:float, dy:float, dz:float):
        """
        Moves the module (dx, dy, dz)
        """
        gmsh.model.occ.translate([ (3,self.id)], dx, dy, dz)