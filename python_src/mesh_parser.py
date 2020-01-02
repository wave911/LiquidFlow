from enum import Enum

class MESH_GEN(Enum):
    SALOME = 0
    GMSH = 1

class SALOME_LINEAR(Enum):
    T1D = 102
    T2D = 203
    T3D = 304
class SALOME_QUAD(Enum):
    T1D = 103
    T2D = 206
    T3D = 310

class MESH_TYPE(Enum):
    T2D = 0
    T3D = 1

class MESH_KEY(Enum):
    Points_Number = "points_number"
    Total_Elements_Number = "total_elements_number"
    Points = "points"
    T1D = "1d"
    T2D = "2d"
    T3D = "3d"

class Point3D(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class Mesh_Parser(object):

    def __init__(self, filepath, mesh_type, mesh_gen=MESH_GEN.SALOME):
        self.filepath = filepath
        self.mesh_type = mesh_type
        self.mesh_gen = mesh_gen
        self.mesh = {}

    def __parse_salome__(self):
        line = None
        try:
            f = open(self.filepath, "r")
            line = f.readline().rstrip().split(' ')
        except Exception as e:
            print("Exception on open {0}. Exception {1}".format(self.filepath, e))
        self.mesh[MESH_KEY.Points_Number] = int(line[0])
        self.mesh[MESH_KEY.Total_Elements_Number] = int(line[1])
        self.mesh[MESH_KEY.Points] = []
        self.mesh[MESH_KEY.T1D] = []
        self.mesh[MESH_KEY.T2D] = []
        self.mesh[MESH_KEY.T3D] = []

        f = self.__get_points__(f)
        f = self.__get_elements__(f)
        return self.mesh

    def __get_points__(self, file):
        points_number = self.mesh[MESH_KEY.Points_Number]
        for i in range(0, points_number):
            line = file.readline().rstrip().split(' ')
            # p = Point3D(float(line[1]), float(line[2]), float(line[3]))
            p = [float(line[1]), float(line[2]), float(line[3])]
            self.mesh[MESH_KEY.Points].append(p)
        return file

    def __get_elements__(self, file):
        lines = file.readlines()
        for line in lines:
            line = line.rstrip().split(' ')
            element_type = int(line[1])
            element = [int(x) - 1 for x in line[2:]]

            if element_type in [int(SALOME_LINEAR.T1D.value), int(SALOME_QUAD.T1D.value)]:
                self.mesh[MESH_KEY.T1D].append(element)
            if element_type in [int(SALOME_LINEAR.T2D.value), int(SALOME_QUAD.T2D.value)]:
                self.mesh[MESH_KEY.T2D].append(element)
            if element_type in [int(SALOME_LINEAR.T3D.value), int(SALOME_QUAD.T3D.value)]:
                self.mesh[MESH_KEY.T3D].append(element)

        pass

    def __parse__(self):
        if self.mesh_gen == MESH_GEN.SALOME:
            self.__parse_salome__()

    def get_mesh(self):
        self.__parse__()
        return self.mesh



