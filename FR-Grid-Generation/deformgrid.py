"""
This Python program will generate a perturbed
grid with gauss lobatto legendre points. It will need as parameters
the minimum and maximum x and y values (for 3D z is needed) and it 
will need the number of grid points to use (for GLL nodes, the 
grid points and solution points coincide with each other)
"""


import numpy
import matplotlib.pyplot as plt
import math

#give the number of elements in each coordinate direction
CONST_DIMENSION = 8

#This gives the order p of the solution approximation. If p =2,
# then 3 solution points are needed in each coordinate direction
CONST_P = 1

#The dimensions of the rectangular grid that will be deformed by adding
# the sinusoidal deformation

CONST_NumProcessors = 1

CONST_xMin = -8.0
CONST_xMax = 8.0
CONST_yMin = -8.0
CONST_yMax = 8.0

CONST_dx = (CONST_xMax-CONST_xMin)/(CONST_DIMENSION)
CONST_dy = (CONST_yMax - CONST_yMin)/(CONST_DIMENSION)

CONST_Deform = True

#CONST_MESHFILENAME = "8x8_P2Riemann_4.msh"

meshfileName = str(CONST_DIMENSION) + "x" + \
    str(CONST_DIMENSION) + "_"
meshSolPtsName = str(CONST_DIMENSION) + "x" + \
    str(CONST_DIMENSION)+"_P"+str(CONST_P)

if(CONST_Deform == True):
     meshfileName = meshfileName + "Sine_" + str(CONST_NumProcessors) + ".msh"
     meshSolPtsName = meshSolPtsName + "Sine_" + str(CONST_NumProcessors)+ "SolPts.msh"
else:
     meshfileName = meshfileName + "Rect_" + str(CONST_NumProcessors) + ".msh"
     meshSolPtsName = meshSolPtsName + "Rect_" + str(CONST_NumProcessors) + "SolPts.msh"

CONST_MESHFILENAMEPERIODIC = meshfileName
CONST_SOLUTIONPTFILENAME = meshSolPtsName

CONST_GaussQuadratureRootsAndCoefficients = {2: [[0.5773502692, -0.5773502692], [1.0, 1.0]],
    3: [[0.7745966692, 0.0000000000, -0.7745966692],
        [0.5555555556, 0.8888888889, 0.5555555556]],
    4: [[0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116],
            [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]]}

CONST_GaussLobattoRootsAndCoefficients = {
    2: [[-1,1],[0,0]],
    3: [[-1., 0, 1.], [0.3333333333333333333333, 1.333333333333333333333, 0.3333333333333333333333]],
    4: [[-1, -0.4472135954999579392818, 0.4472135954999579392818, 1],
        [0.1666666666666666666667, 0.833333333333333333333, 0.833333333333333333333, 0.1666666666666666666667]],
    5: [[-1., -0.654654, 0.0, 0.654654, 1.0],[0,0,0,0,0]]
    }

#This class holds the information for a grid point.
class GridPoint(object):
    def __init__(self,x,y):
        self.x = x
        self.y = y
    

    # The getters
    def getX(self):
        return self.x
    def getY(self):
        return self.y
    def getID(self):
        return self.id
    

    #The setters
    def setX(self,x):
        self.x = x

    def setY(self,y):
        self.y = y


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
        
        """
        if((abs(x - 2.22044e-16) <0.000001) or (abs(y - 2.22044e-16) <0.000001)):
            print " - Map Parent"
            print "     - x: ", x
            print "     - y: ", y
            print "     xi: ", xi
            print "     eta: ", eta
            print " xMaxPhysical: ", xMaxPhysical
            print " xMinPhysical: ", xMinPhysical
            print " yMaxPhysical: ", yMaxPhysical
            print " yMinPhysical: ", yMinPhysical
        """
        
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




#Takes as input the physical element matrix and plotx all the
# elements
def plotElements(PhysicalElementMatrix):
    xVector = []
    yVector = []
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            elementObject = PhysicalElementMatrix[i][j]
            #print " "
            #print "partition: " + str(elementObject.getPartitionNumber())
            #print "xc: " + str(elementObject.getGridPointsMatrix()[CONST_P/2][CONST_P/2].getX())
            #print "yc: " + str(elementObject.getGridPointsMatrix()[CONST_P/2][CONST_P/2].getY())
            for gridPointObject in elementObject.gridPointList:
                xVector.append(gridPointObject.getX())
                yVector.append(gridPointObject.getY())

    plt.scatter(xVector,yVector,s=10,c='b')
    plt.grid()
    # create the lines around all the elements
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            elementObject = PhysicalElementMatrix[i][j]
            xVector,yVector = elementObject.getOuterLineVectors()
            plt.plot(xVector,yVector)

    plt.show(block=True)



#The sinusoidal perturbation function that is used to perturb to
# deform the rectangular mesh. It has as parameters a  and da, which
# is the variable and the grids seperation in that direction. On top of
# this are the parameters n, A and Lo

CONST_n = 2.
CONST_A = 0.5
CONST_Lo = 8.
def sinPertrubFunction(a, da):
    return CONST_A*da*math.sin((CONST_n/CONST_Lo)*3.1415926*a)

#the method used for perturbing the grid.
def perturbGrid(PhysicalElementMatrix):
    # for every element, perturb only the x coordinate using the
    # formula from equation (45-47)
    
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            elementObject = PhysicalElementMatrix[i][j]
            
            for gridObject in elementObject.getGridPoints():
                x = gridObject.getX()
                y = gridObject.getY()
                
                if (CONST_Deform):
                    x = x + sinPertrubFunction(y,CONST_dy)
                    y = y + sinPertrubFunction(x,CONST_dx)
                
                gridObject.setX(x)
                gridObject.setY(y)




# the method that is used for creating the periodic Boundary Conditions
# connections. This information will be stored in tuples (in the order in which
# they will be printed) and then stored into the file.
#for each periodic tuple, the format used to hold the data
# will be (indexA, indexB, faceA, faceB). Here, indexA and indexB
# are the index of the element in the PhysicalElementList + 1. The face
# index will be the index of the face on the element. To find this value,
# find the two nodes that make up the face of the element in question. Then
# the index of the face is found based on the placement of the node values
# in the element's node array

def LoadPeriodicBCData(PhysicalElementMatrix, PhysicalElementList, PeriodicBoundaryConditionsList):
    
    """
    A sample grid
         1
       -----
    2 |      | 4
       -----
         3
    """
    #print "Top Bot: "
    #compute the periodic connections between the bottom and top faces of mesh
    # (i.e. side 1 and 3).
    for i in range(CONST_DIMENSION):
        elementTopRow = PhysicalElementMatrix[0][i]
        elementTopRowIndex = PhysicalElementList.index(elementTopRow)
        elementBotRow = PhysicalElementMatrix[CONST_DIMENSION-1][i]
        elementBotRowIndex = PhysicalElementList.index(elementBotRow)
    
        elementTopRowGPMatrix = elementTopRow.getGridPointsMatrix()
        elementBotRowGPMatrix = elementBotRow.getGridPointsMatrix()
        
        #get the top node points for elementTopRow
        elementTopRowGP1 =elementTopRowGPMatrix[0][0]
        elementTopRowGP2 =elementTopRowGPMatrix[0][CONST_P]
    
        #get the bottom node points for elementBottomRow
        elementBotRowGP1 = elementBotRowGPMatrix[CONST_P][0]
        elementBotRowGP2 = elementBotRowGPMatrix[CONST_P][CONST_P]
    
        #get the face indeces for the two elements
        elementTopRowFaceIndex = elementTopRow.getFaceIndex(elementTopRowGP1, \
                                                            elementTopRowGP2)
        elementBotRowFaceIndex = elementBotRow.getFaceIndex(elementBotRowGP1, \
                                                            elementBotRowGP2)
                                                            
        #store the data in a tuple
        dataTuple = (elementTopRowIndex, elementBotRowIndex, \
                     elementTopRowFaceIndex, elementBotRowFaceIndex)
        
        #print dataTuple
        PeriodicBoundaryConditionsList.append(dataTuple)

    #print "Left Right: "
    #compute the periodic connections between the left and right faces of the mesh
    # (i.e. side 2 and 4).
    for i in range(CONST_DIMENSION):
        elementLeftCol = PhysicalElementMatrix[i][0]
        elementLeftColIndex = PhysicalElementList.index(elementLeftCol)
        elementRightCol = PhysicalElementMatrix[i][CONST_DIMENSION-1]
        elementRightColIndex = PhysicalElementList.index(elementRightCol)
            
        elementLeftColGPMatrix = elementLeftCol.getGridPointsMatrix()
        elementRightColGPMatrix = elementRightCol.getGridPointsMatrix()
            
        #get the Left node points for elementLeftCol
        elementLeftColGP1 =elementLeftColGPMatrix[0][0]
        elementLeftColGP2 =elementLeftColGPMatrix[CONST_P][0]
            
        #get the right node points for elementRightCol
        elementRightColGP1 = elementRightColGPMatrix[0][CONST_P]
        elementRightColGP2 = elementRightColGPMatrix[CONST_P][CONST_P]
            
        #get the face indeces for the two elements
        elementLeftColFaceIndex = elementLeftCol.getFaceIndex(elementLeftColGP1, \
                                                            elementLeftColGP2)
        elementRightColFaceIndex = elementRightCol.getFaceIndex(elementRightColGP1, \
                                                            elementRightColGP2)
                                                                
        #store the data in a tuple
        dataTuple = (elementLeftColIndex, elementRightColIndex, \
                     elementLeftColFaceIndex, elementRightColFaceIndex)

    #print dataTuple
        PeriodicBoundaryConditionsList.append(dataTuple)




# go through all the elements and all those that are on a boundary add to
# the list
def LoadBCNodePointsRiemann(PhysicalElementMatrix, BoundaryConditionNodesList):
    
    # Bottom Boundary
    for i in range(CONST_DIMENSION):
        elementObject = PhysicalElementMatrix[CONST_DIMENSION-1][i]
        GPBotLeft = elementObject.getGridPointsMatrix()[CONST_P][0]
        GPBotRight = elementObject.getGridPointsMatrix()[CONST_P][CONST_P]

        if(getIndexNode(GPBotLeft,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPBotLeft)
        
        if(getIndexNode(GPBotRight,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPBotRight)

    # Top Boundary
    for i in range(CONST_DIMENSION):
        elementObject = PhysicalElementMatrix[0][i]
        GPTopLeft = elementObject.getGridPointsMatrix()[0][0]
        GPTopRight = elementObject.getGridPointsMatrix()[0][CONST_P]
        
        if(getIndexNode(GPTopLeft,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPTopLeft)
        
        if(getIndexNode(GPTopRight,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPTopRight)

    # Left Boundary
    for i in range(CONST_DIMENSION):
        elementObject = PhysicalElementMatrix[i][0]
        GPTopLeft = elementObject.getGridPointsMatrix()[0][0]
        GPBotLeft = elementObject.getGridPointsMatrix()[CONST_P][0]
        
        if(getIndexNode(GPTopLeft,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPTopLeft)

        if(getIndexNode(GPBotLeft,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPBotLeft)

    # Right Boundary
    for i in range(CONST_DIMENSION):
        elementObject = PhysicalElementMatrix[i][CONST_DIMENSION-1]
        GPTopRight = elementObject.getGridPointsMatrix()[0][CONST_P]
        GPBotRight = elementObject.getGridPointsMatrix()[CONST_P][CONST_P]
        
        if(getIndexNode(GPTopRight,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPTopRight)

        if(getIndexNode(GPBotRight,BoundaryConditionNodesList) == -1):
            BoundaryConditionNodesList.append(GPBotRight)


#create a file for all the nodes in the
# right order. That is, place the nodes into the file
# with the node points being printed from bottom to up and
# from left to right go through all the elements and put into
# the node list their bottom left node (for
# the top most element of a column put in the
# top left node also)

def LoadMeshNodePoints(PhysicalElementMatrix, MeshNodePointsList):
    for c in range(CONST_DIMENSION):
        for r in range(CONST_DIMENSION-1, -1 ,-1):
            
            elementObject = PhysicalElementMatrix[r][c]
            elementGridPointMatrix = elementObject.getGridPointsMatrix()
            
            #get the grid point at the bottom left
            elementGridPointbotLeft = elementGridPointMatrix[CONST_P][0]
            
            MeshNodePointsList.append(elementGridPointbotLeft)
            
            
            # if this is the top most element in the column, also
            # add the top left element
            if (r == 0):
                MeshNodePointsList.append(elementGridPointMatrix[0][0])

        # For the last column, also add all the right nodes
    for r in range(CONST_DIMENSION-1, -1 ,-1):
        elementObject = PhysicalElementMatrix[r][CONST_DIMENSION-1]
        MeshNodePointsList.append(elementObject.getGridPointsMatrix()[CONST_P][CONST_P])
        if(r == 0):
            MeshNodePointsList.append(elementObject.getGridPointsMatrix()[0][CONST_P])


# The mesh file that will have periodic boundary conditions
def printMeshFilePERIODIC(PhysicalElementList, MeshNodePointsList, PeriodicBoundaryConditionsList):
    fileName = ""
    #fileName = str(CONST_DIMENSION) + "x" + str(CONST_DIMENSION) + "_4"+"_sine.msh"
    
    fileName = CONST_MESHFILENAMEPERIODIC
    file = open(fileName, "w")
    file.write("Number of grid points: \n")
    file.write(str(len(MeshNodePointsList)) + "\n")
    file.write("Number of QUADS: \n")
    
    file.write(str(len(PhysicalElementList)) + "\n")
    file.write("Nodes coordinates: \n")
    #print the node coordinates in the order in which they appear in the file
    for MeshNodePoint in MeshNodePointsList:
        file.write(str(MeshNodePoint.getX()) + " " + str(MeshNodePoint.getY()) + "\n")
    file.write("Connectivity QUAD: \n")

    for elementObject in PhysicalElementList:
        file.write(elementObject.getConnectivityString() + "\n")
    
    # the riemann boundary condition case
    file.write("PERIODIC " + str(len(PeriodicBoundaryConditionsList)) + "\n")
    for dataTuple in PeriodicBoundaryConditionsList:
        dataTupleString = ""
        for i in range(4):
            dataTupleString = dataTupleString + str(dataTuple[i]) + " "
        file.write(dataTupleString + "\n")

    file.close()


# print everything into the mesh file so that the mesh can
# be read by the higher order code. This mesh file uses Riemann
# boundary conditions
def printMeshFileRIEMANN(PhysicalElementList, MeshNodePointsList, BoundaryConditionNodesList):
    fileName = ""
    fileName = str(CONST_DIMENSION) + "x" + str(CONST_DIMENSION) + "_4"+"_sine.msh"

    #overriding the name of the file. Makes it easier
    # to transfer to guillimin
    fileName = CONST_MESHFILENAME
    file = open(fileName, "w")
    file.write("Number of grid points: \n")
    file.write(str(len(MeshNodePointsList)) + "\n")
    file.write("Number of QUADS: \n")
    
    file.write(str(CONST_DIMENSION**2) + "\n")
    file.write("Nodes coordinates: \n")
    #print the node coordinates in the order in which they appear in the file
    for MeshNodePoint in MeshNodePointsList:
        file.write(str(MeshNodePoint.getX()) + " " + str(MeshNodePoint.getY()) + "\n")
    file.write("Connectivity QUAD: \n")

    for elementObject in PhysicalElementList:
        file.write(elementObject.getConnectivityString() + "\n")
    
    # the riemann boundary condition case
    file.write("RIEMANN " + str(len(BoundaryConditionNodesList)) + "\n")
    for node in BoundaryConditionNodesList:
        #remember to add 1 to each index found
        file.write(str(getIndexNode(node, MeshNodePointsList) + 1) + "\n")

    file.close()


def getIndexNode(NodePoint, NodePointList):
    
    Tolerance = 0.0000001
    for node in NodePointList:
        
        
        if ((node.getX() == NodePoint.getX()) and (node.getY()==NodePoint.getY())):
            return NodePointList.index(node)
        
        
        if((abs(NodePoint.getX() - node.getX())<= Tolerance) and \
            (abs(NodePoint.getY() - node.getY()) <= Tolerance)):
            if (((node.getX() == NodePoint.getX()) and \
                 (node.getY()==NodePoint.getY())) == False):
                return NodePointList.index(node)

    return -1


# Go through each element and print all the solution points. Print
# the points for each element in the same order in which the
# connectivity file has the elements ordered. In addition, for each point, first
# print the node coordinates (for comparing points), and then print all the
# solution points in the exact order they must be loaded into the
# 1D dof array in the higher order code.
def printSolutionPoints(PhysicalElementMatrix):
    
    file = open(CONST_SOLUTIONPTFILENAME, "w")
    
    # find the number of solution points in the file so that C
    # knows to allocate an array of that size to hold all the data. Include
    # in the variable the number of sol pts per element and number of vertices
    numPoints = (CONST_DIMENSION**2)*((CONST_P+1)**2) + (CONST_DIMENSION**2)*4
    
    file.write(str(numPoints) + "\n")
    for c in range(CONST_DIMENSION):
        for r in range(CONST_DIMENSION-1, -1, -1):
            elementObject = PhysicalElementMatrix[r][c]
            elementGridPointMatrix = elementObject.getGridPointsMatrix()
            #first, print the 4 vertices for the element: (from bot right, ccw)
            file.write(str(elementGridPointMatrix[CONST_P][CONST_P].getX()) + " " + \
                       str(elementGridPointMatrix[CONST_P][CONST_P].getY()) + "\n")
            file.write(str(elementGridPointMatrix[0][CONST_P].getX()) + " " + \
                       str(elementGridPointMatrix[0][CONST_P].getY()) + "\n")
            file.write(str(elementGridPointMatrix[0][0].getX()) + " " + \
                       str(elementGridPointMatrix[0][0].getY())+"\n")
            file.write(str(elementGridPointMatrix[CONST_P][0].getX()) + " " + \
                       str(elementGridPointMatrix[CONST_P][0].getY())+"\n")
            
            # print the node points now :
            #The way in which the points must be printed in is by going from bottom
            # to the top from the right of the element to the left of the elem
            
            for c1 in range(CONST_P, -1, -1):
                for r1 in range(CONST_P,-1,-1):
                    GP = elementGridPointMatrix[r1][c1]
                    file.write(str(GP.getX()) + " " + str(GP.getY()) + "\n")

    file.close()




#Go through all the elements and use the already generated MeshNodePointsList
# to create the connectivity string for each element
def ComputeConnectivity(PhysicalElementMatrix, MeshNodePointsList):
    
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            
            elementObject = PhysicalElementMatrix[i][j]
            
            # get the points from the Node array
            elementNodeArray = elementObject.getNodeArray()
            
            # Find the index of all the points and add 1 (the mesh file requires
            # the extra addition by 1)
            elementConnectivityString = ""
            
            for elementNode in elementNodeArray:
                index = getIndexNode(elementNode, MeshNodePointsList) + 1
                if(index == 0):
                    print "         index error: " + str(elementNode.getX()) + \
                        " y: " + str \
                        (elementNode.getY())
                elementConnectivityString = elementConnectivityString + \
                    str(index)+ " "
            elementConnectivityString = elementConnectivityString + \
                    str(elementObject.getPartitionNumber())+ " "+ str(0)
            elementObject.setConnectivityString(elementConnectivityString)


#The main method
def main():
    # Create the vertex points for a test point
    # This matrix will be used for creating the rectangular
    # cartesian grid which will then be deformed.

    #Create the initial rectangular grid. The elements will be placed
    # into a big 2D array which will hold all the element objects. These
    # element objects will be placed starting from the bottom left and then
    # proceed to the top right.
    PhysicalElementMatrix = []
    for i in range(CONST_DIMENSION):
        rowArray = []
        for j in range(CONST_DIMENSION):
            rowArray.append(None)
        PhysicalElementMatrix.append(rowArray)


    # i is for the rows and j is for the columns
    for i in range(CONST_DIMENSION):
        for j in range(CONST_DIMENSION):
            xCenter = CONST_xMin + j*CONST_dx + (CONST_dx)/2.
            yCenter = CONST_yMin + i*CONST_dy + (CONST_dy)/2.
            
            #create a 2D array for the vertices of each element
            VerticesPointsMatrix = []
            for l in range(2):
                rowPointMatrix = []
                for m in range(2):
                    rowPointMatrix.append(None)
                VerticesPointsMatrix.append(rowPointMatrix)
    
    
            #create a physical element about the center
            xMin = xCenter - 0.5*CONST_dx
            xMax = xCenter + 0.5*CONST_dx
            yMin = yCenter - 0.5*CONST_dy
            yMax = yCenter + 0.5*CONST_dy
            
            #fill the above points, which make up the vertices of the element
            # into the verticesPointsMatrix. The points need to be placed
            # in the following format into the verticesMatrix:
            # (xMin,yMax) (xMax,yMax)
            # (xMin,yMin) (xMax,yMin)
            for r in range(2):
                for c in range(2):
                    x = xMin + c*CONST_dx
                    y = yMax - r*CONST_dy
                    
                    #create a point tuple
                    point = (x,y)
                    VerticesPointsMatrix[r][c] = point
        
        
        
            #create the element object using the newly created
            # vertices matrix
            elementObject = element(VerticesPointsMatrix)
                
            #depending on where the element object is on the grid, set its
            # MPI number for which processor will take care of it (this is split
            # among 4 processors so the MPI numbers range from 0 to 3)
            elementRow = CONST_DIMENSION -1 -i
            elementCol = j
           
            if(CONST_NumProcessors == 4): 
           
                elementRowPlus1 = elementRow+1
                elementColPlus1 = elementCol+1
                #now, use which location the element is in to set its partition
                #number
                elementObject.setPartitionNumber(-1)
                if((elementRowPlus1<=CONST_DIMENSION/2) and
                   (elementColPlus1<=CONST_DIMENSION/2)):
                    elementObject.setPartitionNumber(0)
                if((elementRowPlus1>CONST_DIMENSION/2) and
                   (elementColPlus1<=CONST_DIMENSION/2)):
                    elementObject.setPartitionNumber(1)
                if((elementRowPlus1>CONST_DIMENSION/2) and
                   (elementColPlus1>CONST_DIMENSION/2)):
                    elementObject.setPartitionNumber(2)
                if((elementRowPlus1<=CONST_DIMENSION/2) and
                   (elementColPlus1>CONST_DIMENSION/2)):
                    elementObject.setPartitionNumber(3)
            else:
                elementObject.setPartitionNumber(0)
            #arrange the elements into the PhysicalElementMatrix like so
            # elementTopCorner elementTopRightCorner
            # elementBottomCorner elementBottomRightCorner
            PhysicalElementMatrix[elementRow][elementCol] = elementObject

    perturbGrid(PhysicalElementMatrix)





    # a list of all the physical elements. This will be the order in which
    # everything will be printed into the files. Use the Physical element matrix
    # to arrange all the points into the 1D list.
    PhysicalElementList = []
    for c in range(CONST_DIMENSION):
        for r in range(CONST_DIMENSION-1, -1, -1):
            elementObject = PhysicalElementMatrix[r][c]
            PhysicalElementList.append(elementObject)


    MeshNodePointsList = []
    LoadMeshNodePoints(PhysicalElementMatrix, MeshNodePointsList)
    
        # for the Riemann boundary conditions
    #BoundaryConditionNodesList = []
    #LoadBCNodePointsRiemann(PhysicalElementMatrix, BoundaryConditionNodesList)
    
    
        #for the periodic boundary conditions
    PeriodicBoundaryConditionsList = []
    LoadPeriodicBCData(PhysicalElementMatrix, PhysicalElementList, \
                       PeriodicBoundaryConditionsList)



    ComputeConnectivity(PhysicalElementMatrix, MeshNodePointsList)
    #printMeshFileRIEMANN(PhysicalElementList, MeshNodePointsList, BoundaryConditionNodesList)
    printMeshFilePERIODIC(PhysicalElementList, MeshNodePointsList, PeriodicBoundaryConditionsList)
    printSolutionPoints(PhysicalElementMatrix)
    plotElements(PhysicalElementMatrix)


main()












