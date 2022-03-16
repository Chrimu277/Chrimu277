# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt

import Exercise2 as ex2
import Exercise1 as ex1
from scipy.linalg import fractional_matrix_power

#gaussian network model

# V(q) = ..
def potential(q1,q2,qcut=5,k=1):
    r = np.linalg.norm(np.subtract(q1,q2))
    if r < qcut:
        return k/2 * r**2
    else:
        return 0

# returns the positions and masses contained in xyzm_dna.txt
def getFileData(filename):
    q = []  # coordinate list
    m = []  # mass list
    # file öffnen
    f = open(filename, "r")
    lines = f.read().split("\n")    # reads file, splits at newlines

    for i in range(len(lines) - 1):
        tokens = lines[i].split(",")    # splits at commas in 4 tokens
        vec = [(float)(tokens[0]), (float)(tokens[1]), (float)(tokens[2])]
        mass = (float)(tokens[3])
        q.append(vec)   #append data
        m.append(mass)

    return q,m


# makes the Potential-matrix
def makePotentials(q):
    V_ij = [[0.0 for i in range(len(q))] for j in range(len(q))]    # empty matrix
    for i in range(len(q)):
        for j in range(len(q)):
            V_ij[i][j] = potential(q[i], q[j], 5, 1)            # fill

    return V_ij


# creates the hessian matrix
def makeHessian(V,q):
    Hij = V.copy()
    for i in range(len(Hij)):
        for j in range(len(Hij[0])):
            Hij[i][j] = 2*V[i][j] / (q[i]*q[j])
    return Hij


# returns diagonal elements of Matrix
def getDiagElements(M):
    d = [M[i][i] for i in range(len(M))]
    return d

# returns the smallest and the biggest possible eigenvektor
def getMinMaxLambda(m):
    ri = ex1.gershGorinRadiusArray(m)   #gets radius
    diag = np.diagonal(m)
    maxlambs = diag + np.absolute(ri)   # calculate y+,y-
    minlambs = diag - np.absolute(ri)
    maxL = np.max(maxlambs)             # get the max, min
    minL = np.min(minlambs)
    return maxL,minL

# modified deflation method with calculates the smallest eigenvectors, not biggest
def deflation3(matrix, n):
    B = matrix.copy()
    ews = []
    evs = []
    res = []
    for i in range(n):
        I = np.eye(len(B))
        B_dash = np.linalg.inv(B + I)
        evsi, ewsi, resi = ex2.powerMethod(B_dash,1e5,1e-8)

        B = B - ewsi* np.outer(evsi,evsi.transpose().conjugate())
        ews.append(ewsi)
        evs.append(evsi)
        res.append(resi)

    return evs,ews,res


def sortEigenValEigenVek(q_z, eigenVek):

    # copies values out
    Z = q_z.copy()
    eigenVekNew = eigenVek.copy()

    indices = np.argsort(Z)   # gets indices of sorted array
    Z = Z[indices]      # sets new arrays sorted
    eigenVekNew = eigenVekNew[:, indices]

    return Z,eigenVekNew

# --- START -----------------------------------------------------------------------------------------------------------

def main():

    # daten extrahieren aus file
    qs , ms = getFileData("xyzm_dna.txt")
    length = len(qs)
    qs_abs = [np.linalg.norm(qs[i]) for i in range(length)]


    # --- a) ---- create Hessian matrix
    print("\n\n --- a) ------------------------------------------------------")
    # matrix of potentials,  checking the Vij
    Vij = makePotentials(qs)
    print("Vij is \n",np.asmatrix(Vij))
    ex2.checkMatrix(Vij)            # method from exercise 2, needs to be in same project folder


    # --- b) --- describe properties of H and get range of eigenvectors
    # making and checking Hij   -> it is tri-diagonal
    print("\n\n --- b) ------------------------------------------------------")
    Hij = makeHessian(Vij,qs_abs)
    print("Hij is \n",np.asmatrix(Hij))
    ex2.checkMatrix(Hij)

    # radius of eigenvalues (diagonal elements are 0 so radii are already the limits)
    ri_hij = ex1.gershGorinRadiusArray(Hij)
    print("\n\nGershgorin radii of Hij: \n")
    #for i in range(len(np.diagonal(Hij))):
    #   print("EigenValue ",i,": ",np.diagonal(Hij)[i],"+- ",ri_hij[i])


    # --- c) --- Make stiffness matrix K and determine range of eigenvalues
    print("\n\n --- c) ------------------------------------------------------")
    M = np.diag(ms)                                     # mass matrix, diagonal contains masses of atoms
    print("\n\nM is ", M)

    r = np.dot(fractional_matrix_power(M,0.5), qs)      # mass reduced matrix
    print("r is ",r)

    K = np.dot(fractional_matrix_power(M,-0.5), np.dot(Hij, fractional_matrix_power(M,-0.5)))  #stiffness matrix
    print("\n\nK is ", K)
    ri_K = ex1.gershGorinRadiusArray(K)
    print("\n\nGershgorin radii of K: \n")
    for i in range(len(np.diagonal(K))):
        print("EigenValue ", i, ": ", np.diagonal(K)[i], "+- ", ri_K[i])

    # --- d) ---- construct from K a matrix K' for which : eigenvalues < 1,
    print("\n\n --- d) ------------------------------------------------------")
    K_dash = K.copy()
    maxLambda, minLambda = getMinMaxLambda(K_dash)

    while maxLambda > 1 or minLambda < -1:
        K_dash = K_dash*0.95

    # now all eigenvalues are abs(Lambda) < 1
    print(maxLambda,minLambda," are now max, min of eigenvalues")

    # --- e) ---- show that the smallest eigenvalue can be obtained by applying powermethod to (K'-I)^-1
    print("\n\n --- e) ------------------------------------------------------")
    I = np.eye(len(K_dash))                             # einheitsmatrix
    toPower = np.linalg.inv(K_dash + I)
    EVs, EW, Res = ex2.powerMethod(toPower,1e5,1e-5)
    print("smallest EW: ",EW)

    # --- f) --- prove the following equality
    print("\n\n --- f) ------------------------------------------------------")
    p = 0
    eps = 1e-8
    sum = np.zeros(np.shape(K_dash))
    while np.abs(np.max(np.linalg.matrix_power(K_dash,p))) > eps:
        sum += np.linalg.matrix_power(-K_dash,p).dot(K_dash+I)
        p += 1

    print("sum is ", np.asmatrix(sum))
    print("K' is ",np.asmatrix(K_dash))

    # --- g) --- use power method to obtain 10 smallest eigenvectors
    evs, ews, residua = deflation3(K_dash, 10)  # get n smallest eigenvalues and eigenvectors
    print("\n\n10 smallest eigenvalues: \n",np.asmatrix(ews))
    print("\n\n10 corresponding eigenvectors: \n", np.asmatrix(evs))
    print("\n\nresidua: \n", residua)

    print("sizes: ",len(evs), len(np.arange(len(evs))))

    # --- h) --- plot eigenvectors as function of z, individually
    fig = plt.figure()
    colors = ["aqua","black","red","green","yellow","orange","purple","fuchsia","silver","grey"]
    for i in range(len(evs)):
        plt.plot(np.arange(len(evs[0])), evs[i], colors[i], markersize = 1)
    plt.title("Eigenvectors of K' as function of z")
    plt.show()

    # --- i) --- use powermethod to calculate 10 largest eigenvalues of K'
    evs, ews, residua = deflation3(K, 10)  # get n biggest eigenvalues and eigenvectors
    for i in range(len(evs)):
        plt.plot(np.arange(len(evs[0])), evs[i], colors[i], markersize = 1)
    plt.title("Eigenvectors of K as function of z")
    plt.show()


if __name__ == '__main__':
    main();
