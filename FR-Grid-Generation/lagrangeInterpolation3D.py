
"""
This python script is for testing the super and sub 
parametrization. Given a set of node points, it will compute the 
shape functions. In addition, if a parameter is specified at the 
node points, then the interpolated value of this parameter will be computed
at a specified point.
    
"""



"""
    (x,y,z): (7.000000,7.000000,8.250000) (xi,eta,zeta): (-1.000000,-1.000000,-1.000000)
    det: 0.121811
    met1: -0.062500, met2: -0.062500, met3: -0.250000, met4: 0.012276
    (x,y,z): (8.000000,7.000000,8.000000) (xi,eta,zeta): (-1.000000,-1.000000,1.000000)
    det: 0.126264
    met1: -0.062500, met2: -0.000000, met3: -0.250000, met4: -0.004311
    (x,y,z): (7.000000,8.000000,8.000000) (xi,eta,zeta): (-1.000000,1.000000,-1.000000)
    det: 0.125809
    met1: -0.000000, met2: -0.062500, met3: -0.250000, met4: 0.000000
    (x,y,z): (8.000000,8.000000,8.000000) (xi,eta,zeta): (-1.000000,1.000000,1.000000)
    det: 0.125000
    met1: -0.000000, met2: -0.000000, met3: -0.250000, met4: 0.000000
    (x,y,z): (7.250000,7.196424,7.163909) (xi,eta,zeta): (1.000000,-1.000000,-1.000000)
    det: 0.117139
    met1: -0.077749, met2: -0.073003, met3: -0.246807, met4: 0.083676
    (x,y,z): (8.250000,6.931025,6.927392) (xi,eta,zeta): (1.000000,-1.000000,1.000000)
    det: 0.127013
    met1: -0.058390, met2: 0.003370, met3: -0.250656, met4: 0.067089
    (x,y,z): (7.000000,8.250000,6.931025) (xi,eta,zeta): (1.000000,1.000000,-1.000000)
    det: 0.125691
    met1: 0.003612, met2: -0.053910, met3: -0.247769, met4: 0.062500
    (x,y,z): (8.000000,8.000000,7.000000) (xi,eta,zeta): (1.000000,1.000000,1.000000)
    det: 0.125809
    met1: 0.022971, met2: 0.022463, met3: -0.251619, met4: 0.062500
    
    
    
    
    Interp (x,y,z): (7.000000,7.000000,8.250000) (xi,eta,zeta): (-1.000000,-1.000000,-1.000000)
    Interp det: 0.121811
    Interp met1: -0.062500, met2: -0.062500, met3: -0.250000, met4: 0.012276
    Interp (x,y,z): (7.500000,7.000000,8.135299) (xi,eta,zeta): (-1.000000,-1.000000,0.000000)
    Interp det: 0.124037
    Interp met1: -0.062500, met2: -0.031250, met3: -0.250000, met4: 0.003983
    Interp (x,y,z): (8.000000,7.000000,8.000000) (xi,eta,zeta): (-1.000000,-1.000000,1.000000)
    Interp det: 0.126264
    Interp met1: -0.062500, met2: -0.000000, met3: -0.250000, met4: -0.004311
    Interp (x,y,z): (7.000000,7.500000,8.135299) (xi,eta,zeta): (-1.000000,0.000000,-1.000000)
    Interp det: 0.123810
    Interp met1: -0.031250, met2: -0.062500, met3: -0.250000, met4: 0.006138
    Interp (x,y,z): (7.500000,7.500000,8.073223) (xi,eta,zeta): (-1.000000,0.000000,0.000000)
    Interp det: 0.124721
    Interp met1: -0.031250, met2: -0.031250, met3: -0.250000, met4: 0.001991
    Interp (x,y,z): (8.000000,7.500000,8.000000) (xi,eta,zeta): (-1.000000,0.000000,1.000000)
    Interp det: 0.125632
    Interp met1: -0.031250, met2: -0.000000, met3: -0.250000, met4: -0.002155
    Interp (x,y,z): (7.000000,8.000000,8.000000) (xi,eta,zeta): (-1.000000,1.000000,-1.000000)
    Interp det: 0.125809
    Interp met1: -0.000000, met2: -0.062500, met3: -0.250000, met4: 0.000000
    Interp (x,y,z): (7.500000,8.000000,8.000000) (xi,eta,zeta): (-1.000000,1.000000,0.000000)
    Interp det: 0.125405
    Interp met1: -0.000000, met2: -0.031250, met3: -0.250000, met4: 0.000000
    Interp (x,y,z): (8.000000,8.000000,8.000000) (xi,eta,zeta): (-1.000000,1.000000,1.000000)
    Interp det: 0.125000
    Interp met1: -0.000000, met2: -0.000000, met3: -0.250000, met4: 0.000000
    Interp (x,y,z): (7.135299,7.120186,7.700154) (xi,eta,zeta): (0.000000,-1.000000,-1.000000)
    Interp det: 0.119475
    Interp met1: -0.070124, met2: -0.067752, met3: -0.248403, met4: 0.047976
    Interp (x,y,z): (7.635299,7.054061,7.595561) (xi,eta,zeta): (0.000000,-1.000000,0.000000)
    Interp det: 0.123057
    Interp met1: -0.065285, met2: -0.033033, met3: -0.249366, met4: 0.039683
    Interp (x,y,z): (8.135299,6.979706,7.461908) (xi,eta,zeta): (0.000000,-1.000000,1.000000)
    Interp det: 0.126638
    Interp met1: -0.060445, met2: 0.001685, met3: -0.250328, met4: 0.031389
    Interp (x,y,z): (7.073223,7.627299,7.595988) (xi,eta,zeta): (0.000000,0.000000,-1.000000)
    Interp det: 0.122612
    Interp met1: -0.034159, met2: -0.062978, met3: -0.248644, met4: 0.039613
    Interp (x,y,z): (7.573223,7.562942,7.555356) (xi,eta,zeta): (0.000000,0.000000,0.000000)
    Interp det: 0.124317
    Interp met1: -0.029319, met2: -0.028260, met3: -0.249606, met4: 0.035466
    Interp (x,y,z): (8.073223,7.489002,7.488773) (xi,eta,zeta): (0.000000,0.000000,1.000000)
    Interp det: 0.126021
    Interp met1: -0.024480, met2: 0.006458, met3: -0.250569, met4: 0.031319
    Interp (x,y,z): (7.000000,8.135299,7.462501) (xi,eta,zeta): (0.000000,1.000000,-1.000000)
    Interp det: 0.125750
    Interp met1: 0.001806, met2: -0.058205, met3: -0.248885, met4: 0.031250
    Interp (x,y,z): (7.500000,8.073223,7.489002) (xi,eta,zeta): (0.000000,1.000000,0.000000)
    Interp det: 0.125577
    Interp met1: 0.006646, met2: -0.023487, met3: -0.249847, met4: 0.031250
    Interp (x,y,z): (8.000000,8.000000,7.500000) (xi,eta,zeta): (0.000000,1.000000,1.000000)
    Interp det: 0.125405
    Interp met1: 0.011486, met2: 0.011232, met3: -0.250809, met4: 0.031250
    Interp (x,y,z): (7.250000,7.196424,7.163909) (xi,eta,zeta): (1.000000,-1.000000,-1.000000)
    Interp det: 0.117139
    Interp met1: -0.077749, met2: -0.073003, met3: -0.246807, met4: 0.083676
    Interp (x,y,z): (7.750000,7.068975,7.065139) (xi,eta,zeta): (1.000000,-1.000000,0.000000)
    Interp det: 0.122076
    Interp met1: -0.068070, met2: -0.034817, met3: -0.248731, met4: 0.075382
    Interp (x,y,z): (8.250000,6.931025,6.927392) (xi,eta,zeta): (1.000000,-1.000000,1.000000)
    Interp det: 0.127013
    Interp met1: -0.058390, met2: 0.003370, met3: -0.250656, met4: 0.067089
    Interp (x,y,z): (7.135299,7.722074,7.068011) (xi,eta,zeta): (1.000000,0.000000,-1.000000)
    Interp det: 0.121415
    Interp met1: -0.037068, met2: -0.063457, met3: -0.247288, met4: 0.073088
    Interp (x,y,z): (7.635299,7.599891,7.043666) (xi,eta,zeta): (1.000000,0.000000,0.000000)
    Interp det: 0.123913
    Interp met1: -0.027389, met2: -0.025270, met3: -0.249213, met4: 0.068941
    Interp (x,y,z): (8.135299,7.462501,6.978272) (xi,eta,zeta): (1.000000,0.000000,1.000000)
    Interp det: 0.126411
    Interp met1: -0.017710, met2: 0.012916, met3: -0.251138, met4: 0.064794
    Interp (x,y,z): (7.000000,8.250000,6.931025) (xi,eta,zeta): (1.000000,1.000000,-1.000000)
    Interp det: 0.125691
    Interp met1: 0.003612, met2: -0.053910, met3: -0.247769, met4: 0.062500
    Interp (x,y,z): (7.500000,8.135299,6.979706) (xi,eta,zeta): (1.000000,1.000000,0.000000)
    Interp det: 0.125750
    Interp met1: 0.013292, met2: -0.015723, met3: -0.249694, met4: 0.062500
    Interp (x,y,z): (8.000000,8.000000,7.000000) (xi,eta,zeta): (1.000000,1.000000,1.000000) 
    Interp det: 0.125809 
    Interp met1: 0.022971, met2: 0.022463, met3: -0.251619, met4: 0.062500
    
"""





# The node points on the computational domain
"""
CONST_PointsMatrix = [[(-1.,1.),(1.,1.)],\
                [(-1.,-1.),(1.,-1.)]]
"""

CONST_PointsMatrix = []

# The number of nodes along each coordinate direction for the
# domain we are interpolating from.
CONST_Dim = 2

# The order of the domain we are interpolating to.
CONST_PNew = 3

# The value at each node point to be interpolated

CONST_ValuesMatrix = [[[0.121811,0.126264],[0.125809,0.125000]],\
                      [[0.117139,0.127013],[0.125691,0.125809]]]

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
    for i in range(CONST_Dim):
        for j in range(CONST_Dim):
            for k in range(CONST_Dim):
                
                # For this point in the nodes matrix, compute the
                # constant xi lagrange polynomial first
                
                xiNodes = []
                # change the i values (the xi values) while keeping j and k fixed
                for i1 in range(CONST_Dim):
                    xiNodes.append(CONST_PointsMatrix[i1][j][k][0])
                
                etaNodes = []
                # change the j values (the eta values) while keeping the i and k fixed
                for j1 in range(CONST_Dim):
                    etaNodes.append(CONST_PointsMatrix[i][j1][k][1])
                
                zetaNodes = []
                # change the k values (the zeta values) while keeping i and j fixed
                for k1 in range(CONST_Dim):
                    zetaNodes.append(CONST_PointsMatrix[i][j][k1][2])
                
                
                
                
                # Compute the lagrange bases:
                LXi = lambda xi, i=i: (LagrangeBasis_i(xiNodes, i, xi))
                LEta = lambda eta, j=j: (LagrangeBasis_i(etaNodes, j, eta))
                LZeta = lambda zeta, k=k: (LagrangeBasis_i(zetaNodes, k, zeta))
                
                
                # Compute the shape function at the node:
                M = lambda xi,eta,zeta, LXi=LXi,LEta=LEta,LZeta=LZeta: (LXi(xi)* \
                                                                        LEta(eta)*LZeta(zeta))
        
                # store the shape function at the node location
                ShapeFunctionsMatrix[i][j][k] = M


def interpolate(xi,eta,zeta,ShapeFunctionsMatrix):
    solution = 0.0
    for i in range(CONST_Dim):
        for j in range(CONST_Dim):
            for k in range(CONST_Dim):
                value = CONST_ValuesMatrix[i][j][k]
                solution = solution + value*ShapeFunctionsMatrix[i][j][k](xi,eta,zeta)
    
    return solution


def main():
    for i in range(CONST_Dim):
        planeArray = []
        for j in range(CONST_Dim):
            rowArray = []
            for k in range(CONST_Dim):
                rowArray.append(None)
            planeArray.append(rowArray)
        CONST_PointsMatrix.append(planeArray)
    
    # Create the matrix for the Shape functions
    ShapeFunctionsMatrix = []
    for i in range(CONST_Dim):
        planeArray = []
        for j in range(CONST_Dim):
            rowArray = []
            for k in range(CONST_Dim):
                rowArray.append(None)
            planeArray.append(rowArray)
        ShapeFunctionsMatrix.append(planeArray)


    dx = (2.)/(CONST_Dim-1)
    dy = (2.)/(CONST_Dim-1)
    dz = (2.)/(CONST_Dim-1)
    # Compute the points matrix based on the number of sol pts along a coordinate direction:
    for i in range(CONST_Dim):
        for j in range(CONST_Dim):
            for k in range(CONST_Dim):
                x = -1. + dx*i
                y = -1. + dy*j
                z = -1. + dz*k
                
                point = (x,y,z)
                CONST_PointsMatrix[i][j][k] = point
    
    
    print CONST_PointsMatrix

            
    computeShapeFunctions(ShapeFunctionsMatrix)

    # With the shape functions at the node points, interpolate the
    # values at the node points at the new node points for a different
    # solution order.
    
    while(True):
        i = input("give i: ")
        j = input("give j: ")
        k = input("give k: ")
    
        i = int(i)
        j = int(j)
        k = int(k)
    
        M = ShapeFunctionsMatrix[i][j][k]
    
        xi = input("xi: ")
        eta = input("eta: ")
        zeta = input("zeta: ")
    
        print("interp",M(xi,eta,zeta))
        break
    

    for i in range(CONST_PNew):
        for j in range(CONST_PNew):
            for k in range(CONST_PNew):
                xi = CONST_GaussLobattoRootsAndCoefficients[CONST_PNew][0][i]
                eta = CONST_GaussLobattoRootsAndCoefficients[CONST_PNew][0][j]
                zeta = CONST_GaussLobattoRootsAndCoefficients[CONST_PNew][0][k]
            
                print "xi,eta,zeta: ", (xi,eta,zeta)
                print "     value: ", interpolate(xi,eta,zeta,ShapeFunctionsMatrix)


main()


