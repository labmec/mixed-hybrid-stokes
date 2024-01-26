from dataclasses import dataclass
import json
import gmsh
import sys
import os

@dataclass
class TPZMeshModeling:

    @staticmethod
    def PrintJson(JsonData: dict, FileName: str)->None:
        """ 
        PrintJson(JsonData, FileName)

        Reads the information from the 'JsonData' dictionary and writes them in a Json file
        named 'FileName'.

        Return: None

        Types: 
        - 'JsonData': dict
        - 'FileName': string 
        """

        json_object = json.dumps(JsonData, indent=4)

        with open(FileName+".json", "w") as outfile:
            outfile.write(json_object)

    @staticmethod
    def WriteMeshFiles(FileName: str)->None:
        """
        Creates the .geo and .msh files from a gmsh model. 

        Return: None 

        Types:
        - 'FileName': string
        """

        # gmsh.write(FileName + ".geo_unrolled")
        gmsh.write(FileName + ".msh")

    @staticmethod
    def MoveFiles(fileName: str, JsonNewPath: str, MeshNewPath:str)->None:
        """
        Move the .json, .msh, and .geo_unrolled files 
        to an given directory

        Return: None

        Types:
            - 'fileName': string
            - 'JsonNewPath': string
            - 'MeshNewPath': string
        """
        
        oldJson = fileName + ".json"
        oldmsh = fileName + ".msh"
        oldgeo = fileName + ".geo_unrolled"

        os.rename(oldJson, JsonNewPath+oldJson)
        os.rename(oldmsh, MeshNewPath+oldmsh)
        os.rename(oldgeo, MeshNewPath+oldgeo)

    @staticmethod
    def CreatePoints(PointCoordinates: list[tuple], lc: float, kernel: str = 'occ')->None:
        """
        Creates points from a list of 'PointCoordinates' with elements spaced with lc. 
        Important: currently it only generates points with sequential numbering! 

        Return: None

        Types:
        - PointCoordinates: list of list of float
        - step: float
        """
        
        for coord in PointCoordinates:
            x,y,z = coord
            if kernel == 'occ': 
                gmsh.model.occ.addPoint(x,y,z,lc)
            elif kernel == 'built':
                gmsh.model.geo.addPoint(x,y,z,lc)

    @staticmethod
    def CreateLines(LineIndexes: list[tuple], kernel:str = 'occ')->None:
        """
        Creates lines from a list of 'LineIndexes'. In this list, it must be put 
        the index of each point that belongs to each line
        Important: currently it only generates lines with sequential numbering! 

        Types: 
        - 'LineIndexes': list of list of int
        """

        for index in LineIndexes:
            init, end = index
            if kernel == 'occ':
                gmsh.model.occ.addLine(init, end)
            elif kernel == 'built':
                gmsh.model.geo.addLine(init, end)

    @staticmethod
    def CreateCurveLoops(CurveLoopIndexes: list[tuple], kernel: str = 'occ')->None:
        """
        Creates curve loops from the 'CurveLoopIndex' list. In this list, it must be put the
        index of each line that belongs to each curve loop.

        Important: currently it only creates curve loops with sequential numbering

        Return: None

        Types: 
        - 'CurveLoopIndex': list of int
        """

        for  index in (CurveLoopIndexes):
            if kernel == 'occ':
                gmsh.model.occ.addCurveLoop(index)
            elif kernel == 'built':
                gmsh.model.geo.addCurveLoop(index)

    @staticmethod
    def CreatePlanes(PlaneIndexes: list[list], kernel: str = 'occ')->None:
        """
        Creates plane surfaces from the 'PlaneIndex' list. In this list, it must be put the index 
        of each curve loop that belong to each plane surface. 

        Return: None

        Types:
        - 'PlaneIndexes': list of int
        """
        for index in (PlaneIndexes):
            if kernel == 'occ': 
                gmsh.model.occ.addPlaneSurface(index) 
            elif kernel == 'built':
                gmsh.model.geo.addPlaneSurface(index) 

    @staticmethod
    def CreateSurfaceLoop(SurfaceLoopIndexes: list, kernel: str = 'occ')->None:
        """
        Creates curve loops from the 'SurfaceLoopIndex' list. In this list, it must be put the
        index of each surface that belongs to each surface loop.

        Important: currently it only creates surface loops with sequential numbering

        Return: None

        Types: 
        - 'SurfaceLoopIndex': list of int
        """

        for index in (SurfaceLoopIndexes):
            if kernel == 'occ':
                gmsh.model.occ.addSurfaceLoop(index)
            elif kernel == 'built':
                gmsh.model.geo.addSurfaceLoop(index)

    @staticmethod
    def CreateVolumes(VolumesIndexes: list, kernel: str = 'occ')->None:
        """
        Creates volumes from the 'VolumesIndezes' list. In this list, it must be put the index 
        of each surface loop that belong to each volume. 

        Return: None

        Types:
        - 'VolumesIndexes': list of int
        """

        for index in (VolumesIndexes):
            if kernel == 'occ':
                gmsh.model.occ.addVolume([index])
            elif kernel == 'built':
                gmsh.model.geo.addVolume([index])

    @staticmethod
    def CreatePhysicalGroup(GroupData: list)->None:
        """
        Creates the Physical Group from the 'GroupDimension', 'GroupIndex', 'GroupID', 'GroupName' information.
        This function should use the same information to write the Json file! So that, the FEM simulation will 
        be properly set.

        Return: None

        Types: 
        - 'GroupDimension': list of int
        - 'GroupIndex': list of index
        - 'GroupID': list of int
        - 'GroupName': list of string
        """

        for dimTag, id, BCname in GroupData:
            dimension, tag = dimTag
            gmsh.model.addPhysicalGroup(dimension, tag, tag=id, name=BCname)

    @staticmethod
    def CreateCircles(Xcenter: float, Ycenter: float, Zcenter: float, Radius: float) -> int:
        """
        Creates surface circles from lists of 'Xcenter', 'Ycenter', 'Zcenter', and 'Radius'. 
        It first generates the circle as a drawong element in gmsh. Then, it transforms these circles 
        into curve loops to finally converts them into surfaces. Returns the surface circle TAG

        Return: surfaceCircle

        Types:
        - 'Xcenter': float
        - 'Ycenter': float
        - 'Zcenter': float
        - 'Radius': float
        - 'surfaceCircle': int
        """

        circle = gmsh.model.occ.addCircle(Xcenter,Ycenter,Zcenter,Radius)
        curveLoopCircle = gmsh.model.occ.addCurveLoop([circle])
        surfaceCircle = gmsh.model.occ.addPlaneSurface([curveLoopCircle])
        
        return surfaceCircle

    @staticmethod
    def CreateRectangles(coordinates:tuple[float], sideX:float, sideY:float)->None:
        square_list = []
        for coord in coordinates:
            x, y, z = coord
            square = gmsh.model.occ.addRectangle(x, y, z, sideX, sideY)
            square_list.append(square)

        return square_list

    @staticmethod
    def MakeHoles(object: int, holesList: list, holeDim: int)->None:
        """
        Makes holes in a surface domain. Given a 'domain' tag, the 'holeList' tags, and the 'meshDim', it uses the 
        gmsh module cut to calculate the boolean difference the object domain and the object to be cut from it.

        Return: None

        Types:
        - 'domain': int
        - 'holesList': list of int
        - 'meshDim': int
        """
        
        holesTuple = [(holeDim, hole) for hole in holesList]
        gmsh.model.occ.cut([(holeDim,object)], holesTuple)

    @staticmethod
    def TurnOnLabels(*variables: str):
        """
        Turn on the selected entities' labels. 
            - points
            - curves
            - surfaces
            - volumes
        """

        if 'points' in variables:
            gmsh.option.setNumber("Geometry.Points",1)

        if 'curves' in variables:
            gmsh.option.setNumber("Geometry.Curves",1)
        
        if 'surfaces' in variables:
            gmsh.option.setNumber("Geometry.Surfaces",1)

        if 'volumes' in variables:
            gmsh.option.setNumber("Geometry.Volumes",1)

    def TurnOnNumbering(*variables: str):
        """
        Turn on the selected entities' labels. 
            - points
            - curves
            - surfaces
            - volumes
        """

        if 'points' in variables:
            gmsh.option.setNumber("Geometry.PointNumbers",1)

        if 'curves' in variables:
            gmsh.option.setNumber("Geometry.LineNumbers",1)
        
        if 'surfaces' in variables:
            gmsh.option.setNumber("Geometry.SurfaceNumbers",1)

        if 'volumes' in variables:
            gmsh.option.setNumber("Geometry.VolumeNumbers",1)

    @staticmethod
    def TurnOnNormals(size: int=50)->None:
        """
        Display the normal vectors with 'size'
        """
        gmsh.option.setNumber("Geometry.Normals", size)

    @staticmethod
    def TurnOnTangents(size: int=50)->None:
        """
        Display the tnagent vectors with 'size'
        """
        gmsh.option.setNumber("Geometry.Tangents", size)

    @staticmethod
    def ShowModel(meshDim: int =-1)->None:
        """
        Show the model
        """
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run() 

    @staticmethod
    def SetMeshSize(fieldID, surfacesList, meshSize, bigNumber = 1e12):
        fieldOperator = gmsh.model.mesh.field 

        fieldOperator.add("Constant", fieldID)
        fieldOperator.set_number(fieldID, "IncludeBoundary",1)
        fieldOperator.set_numbers(fieldID, "SurfacesList", surfacesList)
        fieldOperator.set_number(fieldID, "VIn", meshSize)
        fieldOperator.set_number(fieldID, "VOut", bigNumber)

        fieldOperator.setAsBackgroundMesh(fieldID)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)