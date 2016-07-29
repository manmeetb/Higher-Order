

"""
    1
     z,k row
      .. . . . .3
    .         ..
 2. . . . . .4  . 7
  .         .  . y, j row
  .  5      . .
  .         ..
6 . . . . . . 8
 
 i row
x

"""


# This python script will generate the mesh file with the connectivity
# and the solution point file for the 3D grids.

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math


#give the number of elements in each coordinate direction
CONST_DIMENSION_X = 8
CONST_DIMENSION_Y = 8
CONST_DIMENSION_Z = 1

#This gives the order p of the solution approximation. If p =2,
# then 3 solution points are needed in each coordinate direction
CONST_P = 2


#The box dimensions
CONST_xMin = -8.0
CONST_xMax = 8.0
CONST_yMin = -8.0
CONST_yMax = 8.0
CONST_zMin = -0.5
CONST_zMax = 0.5

CONST_dx = (CONST_xMax-CONST_xMin)/(CONST_DIMENSION_X)
CONST_dy = (CONST_yMax - CONST_yMin)/(CONST_DIMENSION_Y)
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
    def getPoint(self):
        #return a point tuple
        return (self.x, self.y, self.z)
    
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
class Element(object):
    def __init__(self, VerticesPointsMatrixTuples):
        # The size of the vertices points matrix is 2x2x2
        self.VerticesPointsMatrixTuples = VerticesPointsMatrixTuples

        self.VerticesPointsMatrixNP = []
        for l in range(2):
            rowArray1 = []
            for l1 in range(2):
                rowArray2 = []
                for l3 in range(2):
                    rowArray2.append(None)
                rowArray1.append(rowArray2)
            self.VerticesPointsMatrixNP.append(rowArray1)
        
        # create grid point objects for all the tuples in the vertices matrix
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    x = self.VerticesPointsMatrixTuples[i][j][k][0]
                    y = self.VerticesPointsMatrixTuples[i][j][k][1]
                    z = self.VerticesPointsMatrixTuples[i][j][k][2]
                    
                    GP = GridPoint(x,y,z)
                    self.VerticesPointsMatrixNP[i][j][k] = GP
        
    
    
        #create a grid point matrix and list that will hold all the grid points
        # objects for the given element
        self.gridPointList= []
        
        # make the 3d matrix array
        self.gridPointMatrix = []
        for i in range(CONST_P + 1):
            rowArray1 = []
            for j in range(CONST_P + 1):
                rowArray2 = []
                for k in range(CONST_P + 1):
                    rowArray2.append(None)
                rowArray1.append(rowArray2)
            self.gridPointMatrix.append(rowArray1)

        
        # The grid point matrix will hold the points in the following form
        # place the tensor product points from the minimum x,y,z location
        # on the cube to the maximum location. So i=0, j=0 and k=0 for the
        # grid point matrix is the point at the min xyz location, and as
        # i is changed (while keeping j and k fixed) we are seeing the
        # grid points along the x axis

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
        #create a 3D tensor product of the 1D GLL points
        
        for i in range(CONST_P+1):
            for j in range(CONST_P+1):
                for k in range(CONST_P+1):
                    xi = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][i]
                    eta = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][j]
                    zeta = CONST_GaussLobattoRootsAndCoefficients[CONST_P+1][0][k]
                
                    (x,y,z) = self.mapParentToPhysicalRectangular(xi,eta,zeta)
                    GridPointObject = GridPoint(x,y,z)
                    self.gridPointList.append(GridPointObject)
                    self.gridPointMatrix[i][j][k] = GridPointObject


    # This method is to find where to put the GLL node points on
    # the physical element.
    
    def mapParentToPhysicalRectangular(self,xi,eta,zeta):
        
        # The min and max xyz values for the element
        xMinPhysical = self.VerticesPointsMatrixTuples[0][0][0][0]
        yMinPhysical = self.VerticesPointsMatrixTuples[0][0][0][1]
        zMinPhysical = self.VerticesPointsMatrixTuples[0][0][0][2]
        
        xMaxPhysical = self.VerticesPointsMatrixTuples[1][1][1][0]
        yMaxPhysical = self.VerticesPointsMatrixTuples[1][1][1][1]
        zMaxPhysical = self.VerticesPointsMatrixTuples[1][1][1][2]
        
        
        deltaxPhysical = xMaxPhysical - xMinPhysical
        deltayPhysical = yMaxPhysical - yMinPhysical
        deltazPhysical = zMaxPhysical - zMinPhysical
        
        # the ratio along each side which tells us where to
        # place the node points
        xFactor = (xi-(-1.))/2.
        yFactor = (eta-(-1.))/2.
        zFactor = (zeta-(-1.))/2.
        
        x = xMinPhysical + (deltaxPhysical)*(xFactor)
        y = yMinPhysical + (deltayPhysical)*(yFactor)
        z = zMinPhysical + (deltazPhysical)*(zFactor)
        
        return (x,y,z)
   
   
   
    # The method that is in charge of placing the node points
    # in the correct order for the connectivity file.
    def setNodeArray(self):
        
        elementObjectMatrix = self.gridPointMatrix
        
        self.NodeArray.append(elementObjectMatrix[0][0][CONST_P]) #1
        self.NodeArray.append(elementObjectMatrix[CONST_P][0][CONST_P]) #2
        self.NodeArray.append(elementObjectMatrix[0][CONST_P][CONST_P]) #3
        self.NodeArray.append(elementObjectMatrix[CONST_P][CONST_P][CONST_P]) #4
        
        
        self.NodeArray.append(elementObjectMatrix[0][0][0]) #5
        self.NodeArray.append(elementObjectMatrix[CONST_P][0][0]) #6
        self.NodeArray.append(elementObjectMatrix[0][CONST_P][0]) #7
        self.NodeArray.append(elementObjectMatrix[CONST_P][CONST_P][0]) #8
    
    
    
    #the setters
    def setPartitionNumber(self,a):
        self.partitionNumber = a
    
    def setConnectivityString(self,s):
        self.ConnectivityString = s
    
    """
    #given the two corner node values, return the
    #index of the face based on the position of the nodes
    #in the node array. The parameters are nodes (grid points)
    def getFaceIndex(self, node1, node2):

        
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
    """
    
    
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
    
    def getVerticesMatrixNP(self):
        return self.VerticesPointsMatrixNP
    
    # the x,y,z coordinates of what makes up the other boundary line
    def getOuterLineVectors(self):
        
        #For holding the collection of 4 curves that will be used for
        # drawing the outline of the element.
        ListXVectors = []
        ListYVectors = []
        ListZVectors = []
        
            # get the line for the i=0 plane
        xVector1 = []
        yVector1 = []
        zVector1 = []

        # bottom edge (i,k min const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[0][l][0]
            xVector1.append(GridPointObject.getX())
            yVector1.append(GridPointObject.getY())
            zVector1.append(GridPointObject.getZ())
        
        #right edge (i min const, j max const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[0][CONST_P][l]
            xVector1.append(GridPointObject.getX())
            yVector1.append(GridPointObject.getY())
            zVector1.append(GridPointObject.getZ())
        
        # top edge (i min const, k max const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[0][CONST_P-l][CONST_P]
            xVector1.append(GridPointObject.getX())
            yVector1.append(GridPointObject.getY())
            zVector1.append(GridPointObject.getZ())
        
        #left edge (i min const, j min const)
        for l in range(CONST_P + 1):
            GridPointObject = self.gridPointMatrix[0][0][CONST_P-l]
            xVector1.append(GridPointObject.getX())
            yVector1.append(GridPointObject.getY())
            zVector1.append(GridPointObject.getZ())

        ListXVectors.append(xVector1)
        ListYVectors.append(yVector1)
        ListZVectors.append(zVector1)
        
            # get the line for the i=max const plane
        xVector2 = []
        yVector2 = []
        zVector2 = []

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][l][0]
            xVector2.append(GridPointObject.getX())
            yVector2.append(GridPointObject.getY())
            zVector2.append(GridPointObject.getZ())

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][CONST_P][l]
            xVector2.append(GridPointObject.getX())
            yVector2.append(GridPointObject.getY())
            zVector2.append(GridPointObject.getZ())

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][CONST_P-l][CONST_P]
            xVector2.append(GridPointObject.getX())
            yVector2.append(GridPointObject.getY())
            zVector2.append(GridPointObject.getZ())

        for l in range(CONST_P + 1):
            GridPointObject = self.gridPointMatrix[CONST_P][0][CONST_P-l]
            xVector2.append(GridPointObject.getX())
            yVector2.append(GridPointObject.getY())
            zVector2.append(GridPointObject.getZ())

        ListXVectors.append(xVector2)
        ListYVectors.append(yVector2)
        ListZVectors.append(zVector2)

            # get the line for the j=0 plane
        xVector3 = []
        yVector3 = []
        zVector3 = []
        
        # bottom edge (j,k min const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[l][0][0]
            xVector3.append(GridPointObject.getX())
            yVector3.append(GridPointObject.getY())
            zVector3.append(GridPointObject.getZ())
        
        #right edge (j min const, i max const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][0][l]
            xVector3.append(GridPointObject.getX())
            yVector3.append(GridPointObject.getY())
            zVector3.append(GridPointObject.getZ())

        # top edge (j min const, k max const)
        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P-l][0][CONST_P]
            xVector3.append(GridPointObject.getX())
            yVector3.append(GridPointObject.getY())
            zVector3.append(GridPointObject.getZ())
        
        #left edge (i min const, j min const)
        for l in range(CONST_P + 1):
            GridPointObject = self.gridPointMatrix[0][0][CONST_P-l]
            xVector3.append(GridPointObject.getX())
            yVector3.append(GridPointObject.getY())
            zVector3.append(GridPointObject.getZ())

        ListXVectors.append(xVector3)
        ListYVectors.append(yVector3)
        ListZVectors.append(zVector3)

            # get the line for the j=max plane
        xVector4 = []
        yVector4 = []
        zVector4 = []

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[l][CONST_P][0]
            xVector4.append(GridPointObject.getX())
            yVector4.append(GridPointObject.getY())
            zVector4.append(GridPointObject.getZ())

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P][CONST_P][l]
            xVector4.append(GridPointObject.getX())
            yVector4.append(GridPointObject.getY())
            zVector4.append(GridPointObject.getZ())

        for l in range(CONST_P+1):
            GridPointObject = self.gridPointMatrix[CONST_P-l][CONST_P][CONST_P]
            xVector4.append(GridPointObject.getX())
            yVector4.append(GridPointObject.getY())
            zVector4.append(GridPointObject.getZ())

        for l in range(CONST_P + 1):
            GridPointObject = self.gridPointMatrix[0][CONST_P][CONST_P-l]
            xVector4.append(GridPointObject.getX())
            yVector4.append(GridPointObject.getY())
            zVector4.append(GridPointObject.getZ())

        ListXVectors.append(xVector4)
        ListYVectors.append(yVector4)
        ListZVectors.append(zVector4)

        return (ListXVectors,ListYVectors, ListZVectors)




#Methods

#The sinusoidal perturbation function that is used to perturb to
# deform the rectangular mesh. It has as parameters a  and da, which
# is the variable and the grids seperation in that direction. On top of
# this are the parameters n, A and Lo

CONST_n = 2.
CONST_A = 0.2
CONST_Lo = 10.
def sinPertrubFunction(a,b,da,db,dc): # cNew = perturbFunction result
    return CONST_A*dc*math.sin((CONST_n/CONST_Lo)*3.1415926*a)*math.sin((CONST_n/CONST_Lo)*3.1415926*b)

#the method used for perturbing the grid.
def perturbGrid(PhysicalElementMatrix):
    # for every element, perturb only the x coordinate using the
    # formula from equation (45-47)
    
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                elementObject = PhysicalElementMatrix[i][j][k]
                
                for gridObject in elementObject.getGridPoints():
                    x = gridObject.getX()
                    y = gridObject.getY()
                    z = gridObject.getZ()
                    
                    #x = x + sinPertrubFunction(y,z,CONST_dy,CONST_dz, CONST_dx)
                    #y = y + sinPertrubFunction(x,z,CONST_dx,CONST_dz,CONST_dy)
                    #z = z + sinPertrubFunction(x,y,CONST_dx,CONST_dy,CONST_dz)
                    
                    gridObject.setX(x)
                    gridObject.setY(y)
                    gridObject.setZ(z)




#Takes as input the physical element matrix and plotx all the
# elements
def plotElements(PhysicalElementMatrix):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Axes3D.plot_wireframe(ax, z, x, y)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    

    xVector = []
    yVector = []
    zVector = []
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                elementObject = PhysicalElementMatrix[i][j][k]
                for gridPointObject in elementObject.gridPointList:
                    xVector.append(gridPointObject.getX())
                    yVector.append(gridPointObject.getY())
                    zVector.append(gridPointObject.getZ())

    #Axes3D.scatter(ax, xVector, yVector, zVector, s=30, c='r')

    #Axes3D.plot_wireframe(ax, zSolid, xSolid, ySolid, rstride = 1, cstride = 1, color="b")



    # create the lines around all the elements
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                elementObject = PhysicalElementMatrix[i][j][k]
                ListXVectors,ListYVectors,ListZVectors = elementObject.getOuterLineVectors()
                for l in range(len(ListXVectors)):
                    plt.plot(ListXVectors[l],ListYVectors[l], ListZVectors[l], c = 'b')

    Axes3D.set_ylim(ax, [-10,10])
    Axes3D.set_xlim(ax, [-10,10])
    Axes3D.set_zlim(ax, [-10, 10])
    

    plt.show(block=True)


#The method for loading all the mesh's node points in the correct order
# into the list.
def LoadMeshNodePoints(PhysicalElementMatrix, MeshNodePointsList):
    
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                # for this fixed column in the 3D mesh,
                # add the nodes that are at the min x and y locations
                # from the min z to max z value (i1=0, j1=0, k1 changes)

                elementObject = PhysicalElementMatrix[i][j][k]
                elementObjectMatrix = elementObject.getVerticesMatrixNP()
                
                for k1 in range(2):
                    MeshNodePointsList.append(elementObjectMatrix[0][0][k1])
                
            if(j == CONST_DIMENSION_Y-1):
                # looking down from the top, if we are at the top row in the
                # y direction, then also add the other row of grid points

                elementObject = PhysicalElementMatrix[i][j][k]
                elementObjectMatrix = elementObject.getVerticesMatrixNP()
                
                for k1 in range(2):
                    MeshNodePointsList.append(elementObjectMatrix[0][1][k1])

        if(i == CONST_DIMENSION_X-1):
            # if we are at the last row in the x direction, then we need to add the points
            # on the opposite side of the element
            
            for j in range(CONST_DIMENSION_Y):
                for k in range(CONST_DIMENSION_Z):
                    elementObject = PhysicalElementMatrix[i][j][k]
                    elementObjectMatrix = elementObject.getVerticesMatrixNP()
                    
                    for k1 in range(2):
                        MeshNodePointsList.append(elementObjectMatrix[1][0][k1])

                if(j==CONST_DIMENSION_Y-1):
                    #at the top row in the y direction
                    elementObject = PhysicalElementMatrix[i][j][k]
                    elementObjectMatrix = elementObject.getVerticesMatrixNP()
                
                    for k1 in range(2):
                        MeshNodePointsList.append(elementObjectMatrix[1][1][k1])


# search for the index of the node in a list that holds grid point objects
def getIndexNodePointInNodePointsList(NodePoint, NodePointsList):
    for node in NodePointsList:
        if ((node.getX() == NodePoint.getX()) and (node.getY()==NodePoint.getY()) and \
            (node.getZ() == NodePoint.getZ())):
            return NodePointsList.index(node)

    return -1


# This method is for adding the elements into a list in the exact order in which
# they will be printed into the connectivity file. Print the elements in the same order
# in which the node points are stored in their 1D list
def LoadPhysicalElementList(PhysicalElementMatrix, PhysicalElementList):
    
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                PhysicalElementList.append(PhysicalElementMatrix[i][j][k])


#Go through all the elements and use the already generated MeshNodePointsTuplesList
# to create the connectivity string for each element
def ComputeConnectivity(PhysicalElementMatrix, MeshNodePointsList):
    
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                elementObject = PhysicalElementMatrix[i][j][k]
            
                # get the points from the Node array
                elementNodeArray = elementObject.getNodeArray()
                
                # Find the index of all the points and add 1 (the mesh file requires
                # the extra addition by 1)
                elementConnectivityString = ""
                
                for elementNode in elementNodeArray:
                    index = getIndexNodePointInNodePointsList(elementNode, \
                                                              MeshNodePointsList) + 1
                    elementConnectivityString = elementConnectivityString + \
                        str(index)+ " "
                elementConnectivityString = elementConnectivityString + \
                    str(elementObject.getPartitionNumber())+ " "+ str(0)
                elementObject.setConnectivityString(elementConnectivityString)




# go through all the elements on the boundary and those nodes that are on the
# boundary are added to the list
def LoadBCNodePointsRiemann(PhysicalElementMatrix, BoundaryConditionNodesList):
    
    # i=0 plane
    for j in range(CONST_DIMENSION_Y):
        for k in range(CONST_DIMENSION_Z):
            
            elementObject = PhysicalElementMatrix[0][j][k]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the i1=0 plane
            # on the element to the BCNodesList if it hasn't been
            # added yet
            
            for j1 in range(2):
                for k1 in range(2):
                    NP = elementVerticesMatrix[0][j1][k1]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN i0 plane: " + str(NP.getPoint())
                        BoundaryConditionNodesList.append(NP)

    # i=Max plane
    for j in range(CONST_DIMENSION_Y):
        for k in range(CONST_DIMENSION_Z):
            elementObject = PhysicalElementMatrix[CONST_DIMENSION_X-1][j][k]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the i1=max plane for the element
            
            for j1 in range(2):
                for k1 in range(2):
                    NP = elementVerticesMatrix[1][j1][k1]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN imax plane: " + str(NP.getPoint())
                        BoundaryConditionNodesList.append(NP)

    # j=0 plane
    for i in range(CONST_DIMENSION_X):
        for k in range(CONST_DIMENSION_Z):
            elementObject = PhysicalElementMatrix[i][0][k]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the j1=0 plane for the element
            
            for i1 in range(2):
                for k1 in range(2):
                    NP = elementVerticesMatrix[i1][0][k1]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN j0 plane: " + str(NP.getPoint())
                        BoundaryConditionNodesList.append(NP)

    # j=max plane
    for i in range(CONST_DIMENSION_X):
        for k in range(CONST_DIMENSION_Z):
            elementObject = PhysicalElementMatrix[i][CONST_DIMENSION_Y-1][k]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the j1=max plane for the element
            
            for i1 in range(2):
                for k1 in range(2):
                    NP = elementVerticesMatrix[i1][1][k1]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN jmax plane: " + str(NP.getPoint())
                        BoundaryConditionNodesList.append(NP)


    # k=0 plane
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            elementObject = PhysicalElementMatrix[i][j][0]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the k1=0 plane for the element
            
            for i1 in range(2):
                for j1 in range(2):
                    NP = elementVerticesMatrix[i1][j1][0]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN k0 plane: " + str(NP.getPoint())
                        #BoundaryConditionNodesList.append(NP)

    # k=max plane
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            elementObject = PhysicalElementMatrix[i][j][CONST_DIMENSION_Z-1]
            elementVerticesMatrix = elementObject.getVerticesMatrixNP()
            # add all the node points that are on the k1=max plane for the element
            
            for i1 in range(2):
                for j1 in range(2):
                    NP = elementVerticesMatrix[i1][j1][1]
                    if (getIndexNodePointInNodePointsList(NP, BoundaryConditionNodesList) == -1):
                        print "RIEMANN kmax plane: " + str(NP.getPoint())
                        #BoundaryConditionNodesList.append(NP)




# the method that is used for loading the periodic data for the mesh into tuples.
def LoadPeriodicBCData(PhysicalElementMatrix, PhysicalElementList, PeriodicBCList):
    
    1



def main():
    #make the matrix that will hold all the physical element data. The
    # elements will be placed from the xmin, ymin, zmin location to the
    # xmax,ymax,zmax location
    
    PhysicalElementMatrix = []
    for i in range(CONST_DIMENSION_X):
        rowArray = []
        for j in range(CONST_DIMENSION_Y):
            Array2 = []
            for k in range(CONST_DIMENSION_Z):
                Array2.append(0)
            rowArray.append(Array2)
        PhysicalElementMatrix.append(rowArray)

    #print PhysicalElementMatrix

    # find each grid element
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                xCenter = CONST_xMin + i*CONST_dx + (CONST_dx)/2.
                yCenter = CONST_yMin + j*CONST_dy + (CONST_dy)/2.
                zCenter = CONST_zMin + k*CONST_dz + (CONST_dz)/2.
                
                #print xCenter
                #print yCenter
                #print zCenter
                
                #create a 3D array for the vertices of the element
                VerticesPointsMatrix = []
                for l in range(2):
                    rowArray1 = []
                    for l1 in range(2):
                        rowArray2 = []
                        for l3 in range(2):
                            rowArray2.append(None)
                        rowArray1.append(rowArray2)
                    VerticesPointsMatrix.append(rowArray1)
            
                
                #create a physical element about the center
                xMin = xCenter - 0.5*CONST_dx
                xMax = xCenter + 0.5*CONST_dx
                yMin = yCenter - 0.5*CONST_dy
                yMax = yCenter + 0.5*CONST_dy
                zMin = zCenter - 0.5*CONST_dz
                zMax = zCenter + 0.5*CONST_dz
                
                
                #fill the above points, which make up the vertices of the element
                # into the verticesPointsMatrix. i,j,k=0
                # corresponds to the min xyz and i,j,k = max corresponds to the max
                # xyz
                
                for i1 in range(2):
                    for j1 in range(2):
                        for k1 in range(2):
                            x = xMin + i1*CONST_dx
                            y = yMin + j1*CONST_dy
                            z = zMin + k1*CONST_dz
                            
                            VerticesPointsMatrix[i1][j1][k1] = (x,y,z)
                
                
                #create the element object for this vertices matrix
                elementObject = Element(VerticesPointsMatrix)
                
                #set the element's partition number
                elementObject.setPartitionNumber(0)
                
                PhysicalElementMatrix[i][j][k] = elementObject


    #create the list that will hold all the grid points
    MeshNodePointsList = []
    LoadMeshNodePoints(PhysicalElementMatrix, MeshNodePointsList)

    # add the element objects into a list in the same way/order in which
    # the node points are stored in the node points tuples list
    PhysicalElementList = []
    LoadPhysicalElementList(PhysicalElementMatrix, PhysicalElementList)
    
    for NP in MeshNodePointsList:
        print NP.getPoint()

    ComputeConnectivity(PhysicalElementMatrix, MeshNodePointsList)

    #print the connectivity strings
    for element in PhysicalElementList:
        print element.getConnectivityString()


    # The boundary conditions
    #RiemannBCNodesList = []
    #LoadBCNodePointsRiemann(PhysicalElementMatrix, RiemannBCNodesList)

    PeriodicBCList = []
    LoadPeriodicBCData(PhysicalElementMatrix, PhysicalElementList, PeriodicBCList)



    #perturbGrid(PhysicalElementMatrix)
    #plotElements(PhysicalElementMatrix)



main()







