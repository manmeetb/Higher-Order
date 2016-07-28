
# This python script will generate the mesh file with the connectivity
# and the solution point file for the 3D grids.

import numpy
import matplotlib.pyplot as plt
import math


#give the number of elements in each coordinate direction
CONST_DIMENSION = 2
CONST_DIMENSION_Z = 1

#This gives the order p of the solution approximation. If p =2,
# then 3 solution points are needed in each coordinate direction
CONST_P = 3


#The box dimensions
CONST_xMin = -8.0
CONST_xMax = 8.0
CONST_yMin = -8.0
CONST_yMax = 8.0
CONST_zMin = -8.0
CONST_zMax = 8.0

CONST_dx = (CONST_xMax-CONST_xMin)/(CONST_DIMENSION)
CONST_dy = (CONST_yMax - CONST_yMin)/(CONST_DIMENSION)
CONST_dz = (CONST_zMax-CONST_zMin)/(CONST_DIMENSION_Z)


CONST_GaussQuadratureRootsAndCoefficients = {2: [[0.5773502692, -0.5773502692], [1.0, 1.0]],
    3: [[0.7745966692, 0.0000000000, -0.7745966692],
        [0.5555555556, 0.8888888889, 0.5555555556]],
        4: [[0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116],
            [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]]}

CONST_GaussLobattoRootsAndCoefficients = {
    2: [[-1,1],[0,0]],
    3: [[-1, 0, 1], [0.3333333333333333333333, 1.333333333333333333333, 0.3333333333333333333333]],
    4: [[-1, -0.4472135954999579392818, 0.4472135954999579392818, 1],
        [0.1666666666666666666667, 0.833333333333333333333, 0.833333333333333333333, 0.1666666666666666666667]]}




#The Classes

#This class holds the information for a grid point.
class GridPoint(object):
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    
    # The getters
    def getX(self):
        return self.x
    def getY(self):
        return self.y
    def getZ(self):
        return self.z
    
    #The setters
    def setX(self,x):
        self.x = x
    
    def setY(self,y):
        self.y = y

    def setZ(self,z):
        self.z = z


#This is the data structure for holding all the information about an element.
#It will take as input a vertex which gives the vertices of the element. These vertices
# will then scale the parent element GLL nodes to this physical element.
class element(object):
    def __init__(self, VerticesPointsMatrix):
        # The size of the vertices points matrix is 2x2
        self.VerticesPointsMatrix = VerticesPointsMatrix
        
        #create a grid point matrix and list that will hold all the grid points
        # objects for the given element
        self.gridPointList= []
        
        # make the 2d matrix array
        self.gridPointMatrix = []
        for i in range(CONST_P + 1):
            rowArray = []
            for j in range(CONST_P + 1):
                rowArray.append(None)
            self.gridPointMatrix.append(rowArray)
        
        # The grid point matrix will hold the points in the following form:
        # TopLeft, TopRight
        # BotLeft, BotRight
        self.createGridPoints()
        self.partitionNumber = None
        self.ConnectivityString = ""
        
        #This array will hold the array of grid points for the Node array
        # for the element. So, just like in the FR code, NodeArray[0]
        # will be the same as Triangle -> Node[0]. The ordering has the be consistent
        # with the C FR code
        self.NodeArray = []
        
        self.setNodeArray()
            

    # the method for creating the grid points for the element based
    # on the required order of the solution polynomial
    def createGridPoints(self):
        #create a 2D tensor product of the 1D GLL points
        #print "create grid point"
        for i in range(CONST_P+1):
            for j in range(CONST_P+1):
                xi = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][i]
                eta = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][j]
                
                (x,y) = self.mapParentToPhysicalRectangular(xi,eta)
                GridPointObject = GridPoint(x,y)
                self.gridPointList.append(GridPointObject)
                self.gridPointMatrix[CONST_P-j][i] = GridPointObject
                """
                    print "     i: " + str(CONST_P-j)
                    print "     j: " + str(i)
                    print "     x: " + str(x)
                    print "     y: " + str(y)
                """

    # This method is for mapping from the parent element
    # to a rectangular physical element.
    def mapParentToPhysicalRectangular(self,xi,eta):
        xMinPhysical = self.VerticesPointsMatrix[1][0][0]
        yMinPhysical = self.VerticesPointsMatrix[1][0][1]
        
        #print "xMinPhys: " + str(xMinPhysical)
        #print "yMinPhys: " + str(yMinPhysical)
        
        xMaxPhysical = self.VerticesPointsMatrix[0][1][0]
        yMaxPhysical = self.VerticesPointsMatrix[0][1][1]
        
        #print "xMaxPhys: " + str(xMaxPhysical)
        #print "yMaxPhys: " + str(yMaxPhysical)
        
        deltaxPhysical = xMaxPhysical - xMinPhysical
        deltayPhysical = yMaxPhysical - yMinPhysical
        
        xFactor = (xi-(-1.))/2.
        yFactor = (eta-(-1.))/2.
        
        x = xMinPhysical + (deltaxPhysical)*(xFactor)
        y = yMinPhysical + (deltayPhysical)*(yFactor)
        
        return (x,y)
    
    # The method that is in charge of placing the node points
    # in the correct order. For now, the ordering will be to
    # place the grid points from bottom right and ccw.
    def setNodeArray(self):
        elementObjectMatrix = self.gridPointMatrix
        
        GPBotRight = elementObjectMatrix[CONST_P][CONST_P]
        GPTopRight = elementObjectMatrix[0][CONST_P]
        GPTopLeft = elementObjectMatrix[0][0]
        GPBotLeft = elementObjectMatrix[CONST_P][0]
        
        #load the objects into the NArray
        self.NodeArray.append(GPBotRight)
        self.NodeArray.append(GPTopRight)
        self.NodeArray.append(GPTopLeft)
        self.NodeArray.append(GPBotLeft)
    
    
    #the setters
    def setPartitionNumber(self,a):
        self.partitionNumber = a
    
    def setConnectivityString(self,s):
        self.ConnectivityString = s
    
    #given the two corner node values, return the
    #index of the face based on the position of the nodes
    #in the node array. The parameters are nodes (grid points)
    def getFaceIndex(self, node1, node2):
        
        """
            f=1
            n2 ---- n1
            f=2 |    | f=0
            n3 ---- n0
            f=3
            """
        
        faceIndex = -1
        
        #search through the node array. Break when the first match is found
        # with one of the nodes
        index1 = getIndexNode(node1, self.NodeArray)
        index2 = getIndexNode(node2, self.NodeArray)
        
        index1 = self.NodeArray.index(node1)
        index2 = self.NodeArray.index(node2)
        
        if (index1 <= index2):
            faceIndex = index1
        else:
            faceIndex = index2
        
        #because the nodes loop, add a condition that if the two points
        # are the last and first of the array, then faceIndex shouldn't be 0,
        # it should be the index of the last point of the array
        
        if((index1 == 0 and index2 == (len(self.NodeArray)-1)) or \
           (index2 == 0 and index1 == (len(self.NodeArray)-1))):
            faceIndex = len(self.NodeArray)-1
           
        return faceIndex
    
    
    
    #the getters
    def getNodeArray(self):
        return self.NodeArray
    
    def getConnectivityString(self):
        return self.ConnectivityString
    
    def getPartitionNumber(self):
        return self.partitionNumber
    
    def getGridPoints(self):
        return self.gridPointList
    
    def getGridPointsMatrix(self):
        return self.gridPointMatrix
    
    # the x,y coordinates of what makes up the other boundary line
    def getOuterLineVectors(self):
        xVector = []
        yVector = []
        
        #start from the top left and go down
        
        #left edge
        for i in range(CONST_P + 1):
            GridPointObject = self.gridPointMatrix[i][0]
            xVector.append(GridPointObject.getX())
            yVector.append(GridPointObject.getY())
        
        # bottom edge
        for i in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][i]
            xVector.append(GridPointObject.getX())
            yVector.append(GridPointObject.getY())
        
        # right edge
        for i in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P-i][CONST_P]
            xVector.append(GridPointObject.getX())
            yVector.append(GridPointObject.getY())
        
        # top edge
        for i in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[0][CONST_P-i]
            xVector.append(GridPointObject.getX())
            yVector.append(GridPointObject.getY())
        
        return (xVector,yVector)





#Methods

def main():
    #make the matrix that will hold all the physical element data. The
    # elements will be placed from the xmin, ymin, zmin location to the
    # xmax,ymax,zmax location
    
    PhysicalElementMatrix = []
    for i in range(CONST_DIMENSION):
        rowArray = []
        for j in range(CONST_DIMENSION):
            Array2 = []
            for k in range(CONST_DIMENSION_Z):
                Array2.append(0)
            rowArray.append(Array2)
        PhysicalElementMatrix.append(rowArray)

    print PhysicalElementMatrix

    # find each grid element
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            for k in range(CONST_DIMENSION_Z):
                xCenter = CONST_xMin + i*CONST_dx + (CONST_dx)/2.
                yCenter = CONST_yMin + j*CONST_dy + (CONST_dy)/2.
                zCenter = CONST_zMin + k*CONST_dz + (CONST_dz)/2.
                
                
                
                1

    1


main()







