
"""
This python script is for testing the super and sub 
parametrization. Given a set of node points, it will compute the 
shape functions. In addition, if a parameter is specified at the 
node points, then the interpolated value of this parameter will be computed
at a specified point.
    
"""



# The node points on the computational domain
"""
CONST_PointsMatrix = [[(-1.,1.),(1.,1.)],\
                [(-1.,-1.),(1.,-1.)]]
"""

CONST_PointsMatrix = [[(-1.,1.),(0.,1.),(1.,1.)],\
                      [(-1.,0.),(0.,0.),(1.,0)],\
                      [(-1.,-1.),(0,-1.),(1.,-1.)]]

# The number of nodes along each coordinate direction for the
# domain we are interpolating from.
CONST_Dim = 3

# The order of the domain we are interpolating to.
CONST_PNew = 2

# The value at each node point to be interpolated
"""
CONST_ValuesMatrix = [[0.792893, 0.896447],\
                [1.000000, 1.103553]]
"""
CONST_ValuesMatrix = [[1.232560, 0.941884,1.088508],\
                      [1.070451, 0.934169,1.093762],\
                      [1.134841,1.012240,1.044091]]

CONST_GaussLobattoRootsAndCoefficients = {
    2: [[-1, 1], [0.3333333333333333333333, 0.3333333333333333333333]],
    3: [[-1, 0, 1], [0.3333333333333333333333, 1.333333333333333333333, 0.3333333333333333333333]],
    4: [[-1, -0.4472135954999579392818, 0.4472135954999579392818, 1],
        [0.1666666666666666666667, 0.833333333333333333333, 0.833333333333333333333, 0.1666666666666666666667]]}

# this function computes the lagrange basis at a node
# i. It takes in a list of roots and the index value for
# the root at which the lagrange basis (l_i(x)) is being
# computed
def LagrangeBasis_i(listRoots,index_i,x):
    
    # Compute the numerator of the basis:

    productNumerator = 1.
    # go through each of the roots and take the x value and subtract
    # from it the value of the root. Dont subtract the ith root.
    for i in range(len(listRoots)):
        if (i != index_i):
            productNumerator *= (float(x) - listRoots[i])


    # Compute the denominator of the basis
    productDenominator = 1.
    root_i = listRoots[index_i]
    for i in range(len(listRoots)):
        if (i != index_i):
            productDenominator *= (root_i - listRoots[i])


    solution = float(productNumerator)/float(productDenominator)
    return solution



# The method that is used for computing the shape functions. It will
# then store these functions into a 2D list, where each spot will be the
# function for that corresponding node.
def computeShapeFunctions(ShapeFunctionsMatrix):
    for r in range(CONST_Dim):
        for c in range(CONST_Dim):
            
            # For this point in the nodes matrix, compute the
            # constant xi lagrange polynomial first
            
            xiNodes = []
            # change the cols while keeping row fixed in computational
            # element
            for i in range(CONST_Dim):
                xiNodes.append(CONST_PointsMatrix[r][i][0])
            
            etaNodes = []
            # change the rows while keeping col fixed in computational
            # element
            for j in range(CONST_Dim):
                etaNodes.append(CONST_PointsMatrix[j][c][1])
            
            # Compute the lagrange bases:
            LXi = lambda xi, c=c: (LagrangeBasis_i(xiNodes, c, xi))
            LEta = lambda eta, r=r: (LagrangeBasis_i(etaNodes, r, eta))
            
            # Compute the shape function at the node:
            M = lambda xi,eta, LXi=LXi,LEta = LEta: (LXi(xi)*LEta(eta))
    
            # store the shape function at the node location
            ShapeFunctionsMatrix[r][c] = M


def interpolate(xi,eta,ShapeFunctionsMatrix):
    solution = 0.0
    for i in range(CONST_Dim):
        for j in range(CONST_Dim):
            value = CONST_ValuesMatrix[i][j]
            solution = solution + value*ShapeFunctionsMatrix[i][j](xi,eta)
    
    return solution


def main():
    # Create the matrix for the Shape functions
    ShapeFunctionsMatrix = []
    for i in range(CONST_Dim):
        rowArray = []
        for j in range(CONST_Dim):
            rowArray.append(None)
        ShapeFunctionsMatrix.append(rowArray)
            
    computeShapeFunctions(ShapeFunctionsMatrix)

    # With the shape functions at the node points, interpolate the
    # values at the node points at the new node points for a different
    # solution order.
    for i in range(CONST_PNew):
        for j in range(CONST_PNew):
            xi = CONST_GaussLobattoRootsAndCoefficients[CONST_PNew][0][i]
            eta = CONST_GaussLobattoRootsAndCoefficients[CONST_PNew][0][j]
            
            print "xi,eta: ", (xi,eta)
            print "     value: ", interpolate(xi,eta,ShapeFunctionsMatrix)


main()


