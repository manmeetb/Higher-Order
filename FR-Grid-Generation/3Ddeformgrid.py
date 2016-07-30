

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
CONST_DIMENSION_Z = 8

#This gives the order p of the solution approximation. If p =2,
# then 3 solution points are needed in each coordinate direction
CONST_P = 2


#The output files
CONST_MeshFilenamePeriodic = "8x8x8_P2SinePeriodic_4.msh"
CONST_SolPtsFile = "8x8x8_P2SinePeriodic_4_solpts.msh"

#The box dimensions
CONST_xMin = -8.0
CONST_xMax = 8.0
CONST_yMin = -8.0
CONST_yMax = 8.0
CONST_zMin = -8.0
CONST_zMax = 8.0

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
        
        #print "new element:"
        
        
    
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
        
        # set the matrix that holds the vertices of the element
        self.setVerticesPointsMatrix()
        
        self.partitionNumber = None
        self.ConnectivityString = ""
        
        #This array will hold the array of grid points for the Node array
        # for the element. So, just like in the FR code, NodeArray[0]
        # will be the same as Triangle -> Node[0]. The ordering has the be consistent
        # with the C FR code
        self.NodeArray = []
        self.setNodeArray()

        # The SolPtsList will hodl the solution points in the exact order that they
        # will be loaded into dof[i] arrays for each element.
        self.SolPtsList = []
        self.setSolPtsList()




    # The method in charge of placing the node points in the correct order
    # (the dof[i] order)
    def setSolPtsList(self):
        
        #The ordering of the points is from the xi,eta,zeta = -1,-1,-1 to
        # 1,1,1 (with changing zeta first, then eta then xi). The higher order
        # code can be seen to find what node point on each element corresponds
        # to each xi,eta,zeta value to know the correct ordering.
        
        for k in range(CONST_P,-1,-1):
            for j in range(CONST_P+1):
                for i in range(CONST_P+1):

                    self.SolPtsList.append(self.gridPointMatrix[i][j][k])



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


    # the method for setting the vertices points matrix
    def setVerticesPointsMatrix(self):
        
        GP000 = self.gridPointMatrix[0][0][0]
        GP100 = self.gridPointMatrix[CONST_P][0][0]
        GP110 = self.gridPointMatrix[CONST_P][CONST_P][0]
        GP010 = self.gridPointMatrix[0][CONST_P][0]
        GP001 = self.gridPointMatrix[0][0][CONST_P]
        GP101 = self.gridPointMatrix[CONST_P][0][CONST_P]
        GP111 = self.gridPointMatrix[CONST_P][CONST_P][CONST_P]
        GP011 = self.gridPointMatrix[0][CONST_P][CONST_P]
        
        self.VerticesPointsMatrixNP[0][0][0] = GP000
        self.VerticesPointsMatrixNP[1][0][0] = GP100
        self.VerticesPointsMatrixNP[1][1][0] = GP110
        self.VerticesPointsMatrixNP[0][1][0] = GP010
        self.VerticesPointsMatrixNP[0][0][1] = GP001
        self.VerticesPointsMatrixNP[1][0][1] = GP101
        self.VerticesPointsMatrixNP[1][1][1] = GP111
        self.VerticesPointsMatrixNP[0][1][1] = GP011
        
        
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
    

    #given the corner node values, return the
    #index of the face based on the position of the nodes
    #in the node array. The parameters are nodes in a list(grid points)
    def getFaceIndex(self, nodesList):

        faceIndex = -1
        
        
        
        #search through the node array. Break when the first match is found
        # with one of the nodes
        NodeArrayIndexList = []
        for np in nodesList:
            index = getIndexNodePointInNodePointsList(np, self.NodeArray)
            NodeArrayIndexList.append(index + 1)
        
        # get the face index
        
        if((1 in NodeArrayIndexList) and (2 in NodeArrayIndexList) and \
           (3 in NodeArrayIndexList) and (4 in NodeArrayIndexList)):
           faceIndex = 0
           
        if((5 in NodeArrayIndexList) and (7 in NodeArrayIndexList) and \
              (6 in NodeArrayIndexList) and (8 in NodeArrayIndexList)):
            faceIndex = 1
              
        if((1 in NodeArrayIndexList) and (2 in NodeArrayIndexList) and \
            (5 in NodeArrayIndexList) and (6 in NodeArrayIndexList)):
            faceIndex = 2
            
        if((3 in NodeArrayIndexList) and (4 in NodeArrayIndexList) and \
           (7 in NodeArrayIndexList) and (8 in NodeArrayIndexList)):
            faceIndex = 3

        if((1 in NodeArrayIndexList) and (3 in NodeArrayIndexList) and \
           (5 in NodeArrayIndexList) and (7 in NodeArrayIndexList)):
            faceIndex = 4

        if((2 in NodeArrayIndexList) and (4 in NodeArrayIndexList) and \
           (6 in NodeArrayIndexList) and (8 in NodeArrayIndexList)):
            faceIndex = 5

        return faceIndex
    
    
    #the getters
    def getSolutionPointsList(self):
        return self.SolPtsList
    
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
                    
                    x = x + sinPertrubFunction(y,z,CONST_dy,CONST_dz, CONST_dx)
                    y = y + sinPertrubFunction(x,z,CONST_dx,CONST_dz,CONST_dy)
                    z = z + sinPertrubFunction(x,y,CONST_dx,CONST_dy,CONST_dz)
                    
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
    
    # add the min x and y
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            for k in range(CONST_DIMENSION_Z):
                
                # go along the column
                
                elementObj = PhysicalElementMatrix[i][j][k]
                elementObjVertices = elementObj.getVerticesMatrixNP()
                
                # add only the element at the min x and y position
                np = elementObjVertices[0][0][0]
                MeshNodePointsList.append(np)
                
                # if we are at the top of the column, also add the element on top
                # of np in the z direction
                
                if (k == CONST_DIMENSION_Z-1):
                    np = elementObjVertices[0][0][1]
                    MeshNodePointsList.append(np)
            
            # If we are at last row in the j direction, we need to add the vertices that
            # are on the opposite side in the y direction (the max y, min x and min z)
            
            if (j == CONST_DIMENSION_Y-1):
                for k in range(CONST_DIMENSION_Z):
                    elementObj = PhysicalElementMatrix[i][j][k]
                    elementObjVertices = elementObj.getVerticesMatrixNP()
                    
                    #add the node point
                    np = elementObjVertices[0][1][0]
                    MeshNodePointsList.append(np)
                    
                    #If we are at the top of the column
                    
                    if (k == CONST_DIMENSION_Z-1):
                        np = elementObjVertices[0][1][1]
                        MeshNodePointsList.append(np)
                            
        #If we are at the last row along the i direction
        if (i == CONST_DIMENSION_X-1):
            # add the vertices from the neglected lines thus far
            for j in range(CONST_DIMENSION_Y):
                for k in range(CONST_DIMENSION_Z):
                
                    # go along the column
                
                    elementObj = PhysicalElementMatrix[i][j][k]
                    elementObjVertices = elementObj.getVerticesMatrixNP()
                
                    np = elementObjVertices[1][0][0]
                    MeshNodePointsList.append(np)
                
                    # if we are at the top of the column, also add the element on top
                    # of np in the z direction
                
                    if (k == CONST_DIMENSION_Z-1):
                        np = elementObjVertices[1][0][1]
                        MeshNodePointsList.append(np)
            
                    # If we are at last row in the j direction, we need to add the vertices that
                    # are on the opposite side in the y direction (the max y, min x and min z)
            
                if (j == CONST_DIMENSION_Y-1):
                    for k in range(CONST_DIMENSION_Z):
                        elementObj = PhysicalElementMatrix[i][j][k]
                        elementObjVertices = elementObj.getVerticesMatrixNP()
                    
                        #add the node point
                        np = elementObjVertices[1][1][0]
                        MeshNodePointsList.append(np)
                    
                        #If we are at the top of the column
                    
                        if (k == CONST_DIMENSION_Z-1):
                            np = elementObjVertices[1][1][1]
                            MeshNodePointsList.append(np)




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
                    if (index == 0):
                        print "     Error:"
                        print (elementNode.getPoint())
                    
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
    #the periodic BC tuple has to hold the following information:
    # (index element 1, index element 2, face element 1, face element 2)
    
    
    #compute the periodic BCs between k=0 and k=max plane for the grid
    for i in range(CONST_DIMENSION_X):
        for j in range(CONST_DIMENSION_Y):
            elementk0 = PhysicalElementMatrix[i][j][0]
            elementkMax = PhysicalElementMatrix[i][j][CONST_DIMENSION_Z-1]
            
            #the k1 = max plane for elementkmax and k1 = 0 plane for
            # elementkmin are what are connected
            
            # get the collection of node points that make up the face
            # of question for element k0
            elementk0faceNodeList = []
            for i1 in range(2):
                for j1 in range(2):
                    elementk0faceNodeList.append(elementk0.getVerticesMatrixNP()[i1][j1][0])

            elementk0faceIndex = elementk0.getFaceIndex(elementk0faceNodeList)
                
            # get the collection of node points that make up the face
            # of question for element k0
            elementkMaxfaceNodeList = []
            for i1 in range(2):
                for j1 in range(2):
                    elementkMaxfaceNodeList.append(elementkMax.getVerticesMatrixNP()[i1][j1][1])

            elementkMaxfaceIndex = elementkMax.getFaceIndex(elementkMaxfaceNodeList)

            dataTuple = (PhysicalElementList.index(elementk0), \
                         PhysicalElementList.index(elementkMax), elementk0faceIndex, \
                         elementkMaxfaceIndex)

            PeriodicBCList.append(dataTuple)

    #compute the periodic BCs between the j=0 and j=max plane for the grid
    for i in range(CONST_DIMENSION_X):
        for k in range(CONST_DIMENSION_Z):
            elementj0 = PhysicalElementMatrix[i][0][k]
            elementjMax = PhysicalElementMatrix[i][CONST_DIMENSION_Y-1][k]
            
            
            # get the collection of node points that make up the face
            # of question for element j0
            elementj0faceNodeList = []
            for i1 in range(2):
                for k1 in range(2):
                    elementj0faceNodeList.append(elementj0.getVerticesMatrixNP()[i1][0][k1])
        
            elementj0faceIndex = elementj0.getFaceIndex(elementj0faceNodeList)
            
            # get the collection of node points that make up the face
            # of question for element jmax
            elementjMaxfaceNodeList = []
            for i1 in range(2):
                for k1 in range(2):
                    elementjMaxfaceNodeList.append(elementjMax.getVerticesMatrixNP()[i1][1][k1])
            
            elementjMaxfaceIndex = elementjMax.getFaceIndex(elementjMaxfaceNodeList)
            
            dataTuple = (PhysicalElementList.index(elementj0), \
                         PhysicalElementList.index(elementjMax), elementj0faceIndex, \
                         elementjMaxfaceIndex)
                
            PeriodicBCList.append(dataTuple)


    #compute the periodic BCs between the i=0 and i=max plane for the grid
    for j in range(CONST_DIMENSION_Y):
        for k in range(CONST_DIMENSION_Z):
            elementi0 = PhysicalElementMatrix[0][j][k]
            elementiMax = PhysicalElementMatrix[CONST_DIMENSION_X-1][j][k]
        
        
            # get the collection of node points that make up the face
            # of question for element i0
            elementi0faceNodeList = []
            for j1 in range(2):
                for k1 in range(2):
                    elementi0faceNodeList.append(elementi0.getVerticesMatrixNP()[0][j1][k1])

            elementi0faceIndex = elementi0.getFaceIndex(elementi0faceNodeList)

            # get the collection of node points that make up the face
            # of question for element jmax
            elementiMaxfaceNodeList = []
            for j1 in range(2):
                for k1 in range(2):
                    elementiMaxfaceNodeList.append(elementiMax.getVerticesMatrixNP()[1][j1][k1])
        
            elementiMaxfaceIndex = elementiMax.getFaceIndex(elementiMaxfaceNodeList)
                
            dataTuple = (PhysicalElementList.index(elementi0), \
                        PhysicalElementList.index(elementiMax), elementi0faceIndex, \
                             elementiMaxfaceIndex)
                             
            PeriodicBCList.append(dataTuple)



# The method for printing everything into a file so that it can be loaded into the
# FR code.
def printMeshFilePeriodic(PhysicalElementList, MeshNodePointsList, PeriodicBCList):
    fileName = ""
    #fileName = str(CONST_DIMENSION) + "x" + str(CONST_DIMENSION) + "_4"+"_sine.msh"
    
    fileName = CONST_MeshFilenamePeriodic
    file = open(fileName, "w")
    file.write("Number of grid points: \n")
    file.write(str(len(MeshNodePointsList)) + "\n")
    file.write("Number of HEXAHEDRAL: \n")
    
    file.write(str(len(PhysicalElementList)) + "\n")
    file.write("Nodes coordinates: \n")
    #print the node coordinates in the order in which they appear in the file
    for MeshNodePoint in MeshNodePointsList:
        file.write(str(MeshNodePoint.getX()) + " " + str(MeshNodePoint.getY()) + " " + \
                   str(MeshNodePoint.getZ()) + "\n")

    file.write("Connectivity HEXAHEDRAL: \n")
    
    for elementObject in PhysicalElementList:
        file.write(elementObject.getConnectivityString() + "\n")


    file.write("PERIODIC " + str(len(PeriodicBCList)) + "\n")
    for dataTuple in PeriodicBCList:
        dataTupleString = ""
        for i in range(4):
            dataTupleString = dataTupleString + str(dataTuple[i]) + " "
        file.write(dataTupleString + "\n")

    file.close()



# This method will take in the list of elements that make up
# the mesh and print out the inner solution points. This
# is done by first printing the vertices (so that matches can be found), and then
# printing the solution points that will go into the dof arrays for the elements

def printSolutionPointsFile(PhysicalElementList):
    
    file = open(CONST_SolPtsFile, 'w')
    
    # for each element, there are 8 vertices and CONST_P+1^3 solution points
    numLines = (CONST_DIMENSION_X*CONST_DIMENSION_Y*CONST_DIMENSION_Z)*8 + \
        (CONST_DIMENSION_X*CONST_DIMENSION_Y*CONST_DIMENSION_Z)*(CONST_P+1)**3
    
    file.write(str(numLines) + "\n")
    for elementObj in PhysicalElementList:
        # first print the vertices of the element
        elementObjVertices = elementObj.getVerticesMatrixNP()
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    np = elementObjVertices[i][j][k]
                    file.write(str(np.getX()) + " " + str(np.getY()) + " " + str(np.getZ()) + "\n")

    
        # print out the solution points of the element
        elementObjSolPts = elementObj.getSolutionPointsList()
        for sol in elementObjSolPts:
            file.write(str(sol.getX()) + " " + str(sol.getY()) + " " + str(sol.getZ()) + "\n")

    file.close()



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
                
                #set the element's partition number. We will partition the elements by
                # setting a 2 x 2 grid along the xy plane of the mesh
                
                iPlus1 = i+1    # a number from 1 to 8
                jPlus1 = j+1    # a number from 1 to 8
                
                if((iPlus1<=CONST_DIMENSION_X/2) and # i from 1-4, j from 1-4
                   (jPlus1<=CONST_DIMENSION_Y/2)):
                    elementObject.setPartitionNumber(0)
                if((iPlus1>CONST_DIMENSION_X/2) and # i from 5-8, j from 1-4
                   (jPlus1<=CONST_DIMENSION_Y/2)):
                    elementObject.setPartitionNumber(1)
                if((iPlus1>CONST_DIMENSION_X/2) and # i from 5-8, j from 5-8
                   (jPlus1>CONST_DIMENSION_Y/2)):
                    elementObject.setPartitionNumber(2)
                if((iPlus1<=CONST_DIMENSION_X/2) and # i from 1-4, j from 5-8
                  (jPlus1>CONST_DIMENSION_Y/2)):
                    elementObject.setPartitionNumber(3)
                
                
                #elementObject.setPartitionNumber(0)
                
                PhysicalElementMatrix[i][j][k] = elementObject

    perturbGrid(PhysicalElementMatrix)
    
    #create the list that will hold all the grid points
    MeshNodePointsList = []
    LoadMeshNodePoints(PhysicalElementMatrix, MeshNodePointsList)

    # add the element objects into a list in the same way/order in which
    # the node points are stored in the node points tuples list
    PhysicalElementList = []
    LoadPhysicalElementList(PhysicalElementMatrix, PhysicalElementList)

    # get the connectivity strings for each element
    ComputeConnectivity(PhysicalElementMatrix, MeshNodePointsList)


    # The boundary conditions
    #RiemannBCNodesList = []
    #LoadBCNodePointsRiemann(PhysicalElementMatrix, RiemannBCNodesList)

    PeriodicBCList = []
    LoadPeriodicBCData(PhysicalElementMatrix, PhysicalElementList, PeriodicBCList)


    # print everything into a file
    printMeshFilePeriodic(PhysicalElementList, MeshNodePointsList, PeriodicBCList)
    printSolutionPointsFile(PhysicalElementList)


    plotElements(PhysicalElementMatrix)


main()







