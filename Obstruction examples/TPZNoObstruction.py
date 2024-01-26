"""
Class to create a module without obstructions

Created by Carlos Puga: 01/13/2024
"""

#%% ****************** 
#   IMPORTED MODULES
#   ******************
from dataclasses import dataclass, field
import gmsh

from TPZModuleTypology import TPZModuleTypology
#%% ****************** 
#   CLASS DEFINITION
#   ******************
@dataclass
class TPZNoObstruction(TPZModuleTypology):
    """
    This class creates a module without obstrucitons
    
    It inherits a base creator of typologies (TPZModuleTypology)
    to generates the module body and insert the obstruction.

    Fields: 
        - module_typology: information about the type 
        of module to be constructed. Currently is only availabe rectangular 
        modules. More information are required: 
            * Rectangular -> {dx: length on x-direction, dy: length on y-direction}

        - id: obstruction volume's identification (provided by the class itself)
    """
#   ****************** 
#      INITIALIZOR
#   ******************  
    _module_typology: tuple[str, dict]
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

            self.CreateBox(dx, dy)

        elif m_type == 'Circular':
            if not 'radius' in keys:
                self.DebugStop('ERROR: Not enough information to create the box!')

            radius = m_characteristics['radius']

            self.CreateCylinder(radius)
            
        else:
            self.DebugStop(f"ERROR: I'm sorry the typology {m_type} has not been implemented yet")


    def CreateDomain(self)->None:
        """
        Generates the module with a circular obstruction on gmsh.
        """
        self.CheckTypology()

        domain_surfaces = self.surfaces + [self.obstruction_face]

        # creating the module volume
        surface_loop = gmsh.model.occ.addSurfaceLoop(domain_surfaces)
        domain_volume = gmsh.model.occ.addVolume([surface_loop])

        self.id = domain_volume

    def Move(self, dx:float, dy:float, dz:float)->None:
        """
        Moves the module (dx, dy, dz)
        """
        gmsh.model.occ.translate([ (3,self.id)], dx, dy, dz)