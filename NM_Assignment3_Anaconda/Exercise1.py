# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt
import time


# ----------------------------------------------------------------------------------------------------------------------

# calculate 1 rotation matrix
def U_t(i,j,phi,N):

    U = [[0 for i in range(N)]  for j in range(N)]

    #diagonale auf 1
    for y in range(N):
        for x in range(N):
            if( y == x):
                U[y][x] = 1;

    U[i][i] = np.cos(phi)
    U[j][j] = np.cos(phi)
    U[i][j] = -np.sin(phi)
    U[j][i] = np.sin(phi)
    U[0][0] = 1
    U[N-1][N-1] = 1

    return U


# jacobi's method for symmetrical matrics to find solution of eigenproblem
def jacobi(A,eps = 1e-9,maxCount = 1e9,phifac = 0.01):

    A_new = A.copy() # copy to not edit original, faster
    diagonal = False # flag
    counter = 0     # counts each loop
    U_final = np.eye(len(A),len(A))     #rotation matrix (final)
    n = len(A)
    converges = []
    gersh_rads = []

    while diagonal is False and counter < maxCount:
        counter += 1
        i, j, max = getMaxElementAndIndices(A_new);  # rotiere rund um größtes element !
        phi = phifac * np.arctan(2 * A_new[i,j] / (A_new[i][i] - A_new[j][j]))   # 0.01 instead of 0.5

        U = U_t(i, j, phi, n)
        U_T = np.transpose(U)

        U_final = np.dot(U_final,U)         #U = U0 * U1 * .... Ut  .. enthält eigenvektoren !
        A_new = U_T.dot(np.dot(A_new, U))   #rotierte neue matrix A

        # get gershgorins riadus
        gersh_rads.append(gershGorinRadiusArray(A_new))

        # check convergence
        S = convergenceS(A_new)
        converges.append(S)
        if S < eps:
            diagonal = True;

    #Eigenwerte check, U(t) check ->

    A_new = setNonDiagonalZero(A_new)   # kann man auskommentieren ! bessere lesbarkeit
    return A_new,U_final,counter,converges,gersh_rads;


# creates Rnadom symmetrical matrix
def randomMatrix(N):
    M1 = [[(np.random.normal(0,2)) for i in range(N)]  for j in range(N)]
    return M1+np.transpose(M1)

# creates random complex symmetrical
def randomComplexMatrix(N):
    M1 : complex = [[(np.random.normal(0,2)) + 1j*(np.random.normal(0,2)) for k in range(N)]  for l in range(N)]
    return M1+np.transpose(M1)


# searches trhough matrix for absoulte biggest element and its indices
def getMaxElementAndIndices(A):
    N = len(A)
    Imax = 0;
    Jmax = 0;
    Max = 0;

    # seek for the biggest matrix element (absolute) and get indices
    for i in range(N):
        for j in range(i+1,N):
            val = np.abs(A[i][j])
            if val > Max:
                Max = val
                Imax = i
                Jmax = j
    return [Imax,Jmax,Max]

# get sum of 1 row (i) without diagonal element
def rowAbsSum(i,A):
    sum = 0.0;
    for a in range(len(A)):
        sum += np.abs(A[i][a])
    return sum;

# get sum of 1 col (j) without diagonal element
def colAbsSum(j,A):
    sum = 0.0;
    for a in range(len(A)):
        if a != j:
            sum += np.abs(A[a][j])
    return sum;

# calculate radius of greshgorin disks, return as array
def gershGorinRadiusArray(A):
    n = len(A)
    gersh_rads = np.arange(n, dtype = float)
    for i in range(n):
        colrad = colAbsSum(i,A)     # sum of absolute values of the column
        rowrad = rowAbsSum(i,A)     # ...  ...                       row
        if colrad < rowrad:
            gersh_rads[i] = colrad
        else:
            gersh_rads[i] = rowrad

    return gersh_rads;

# check S if matrix has converged (diagonally), convergence = Sum of all off diagonal elements
def convergenceS(A):
    A_temp = A.copy();
    S = 0.0
    for i in range(len(A_temp)):
        for j in range(len(A_temp)):
            if(i != j): # if off diagonal
                S += (A_temp[i][j])**2      # add the square of the element to the sum
    return S;

# miscellaneous
def setNonDiagonalZero(A):
    A_work = A.copy()
    for i in range(len(A_work)):
        for j in range(len(A_work)):
            if i != j:
                A_work[i][j] = 0;
    return A_work

# sorts eigenvalues ascending, and the corresponding eigenvectors alike
def sortEigenValEigenVek(eigenVal, eigenVek):

    # copies values out
    eigenValNew = eigenVal.copy()
    eigenVekNew = eigenVek.copy()

    indices = np.argsort(eigenValNew)   # gets indices of sorted array
    eigenValNew = eigenValNew[indices]      # sets new arrays sorted
    eigenVekNew = eigenVekNew[:, indices]

    return eigenValNew,eigenVekNew

# ---------------------------------------------------------------------------------------------------------------------
def main():

    # symmetric matrices (many hamiltonians can be expressed so)
    # ---a)--- Write function that generates a random NxN matrix

    print("\n\n--- a) ---------------------")
    random1 = randomMatrix(10)
    print("\n\nrandom matrix size 10x10 =\n",np.asmatrix(random1))


    # ---b)--- write function jacobi, that performs repeated unitary transformations of A
    # until: change of matrix elements is very little, size 5x5, result is diagonal matrix

    _N = 5
    print("\n\n--- b) ---------------------")
    A = randomMatrix(_N)                            # make random matrix
    print("\n\na random matrix A =\n",np.asmatrix(A))

    # erzeuge eigenwerte und eigenvektoren, sortiert
    D,U,count,S,gershrad = jacobi(A,1e-9,1e5,0.5)
    eigenvals = np.arange(len(A),dtype=float)
    print("\n\nrotated with jacobi : A' bzw. D =\n",D)
    print("\n\ntransformation matrix U =\n",U)

    gershrad_new = np.transpose(gershrad)
    # plot gershgorin radii
    plt.figure()
    legends = [("radius "+str(i)) for i in range(len(gershrad[0]))]
    plt.plot(np.arange(len(gershrad)),gershrad_new[0])
    plt.plot(np.arange(len(gershrad)),gershrad_new[1])
    plt.plot(np.arange(len(gershrad)),gershrad_new[2])
    plt.plot(np.arange(len(gershrad)),gershrad_new[3])
    plt.plot(np.arange(len(gershrad)),gershrad_new[4])
    plt.title("Gershgorin radii along the iteration steps")
    plt.xlabel("steps")
    plt.ylabel("ri")
    plt.legend(legends)
    plt.show()


    # print counts
    print("\n\nsteps needed:", count)
    # jacobi rotation successfull

    # get eigenvalues
    for i in range(len(D)):
        eigenvals[i] = D[i][i]

    #sort
    eigenvals,U = sortEigenValEigenVek(eigenvals,U)

    #print eigenvals
    print("\n\neigenvalues (found) : \n",eigenvals)
    # print eigenvektors
    print("eigenvektors (found) : \n")
    for i in range(len(U)):
        print("e", i, " = ", U[i])


    # ---c)--- calculate eigenvectors and eigenvalues of A------------------------------------------------------
    print("\n\n--- c) ---------------------")
    eigen = np.linalg.eig(A)
    #sort
    eigenValSort,eigenVekSort = sortEigenValEigenVek(eigen[0],eigen[1])

    print("\n\neigenvalues (builtin method) : \n",eigenValSort)
    print("eigenvektors (builtin method) : \n")
    for i in range(len(eigenVekSort)):
        print("e",i," = ",eigenVekSort[i])


    # ---d)--- apply gershgorin disc method to random symmetric matrix 5x5-----------------------------------------

    print("\n\n--- d) ---------------------")
    _N = 5
    random2 = randomComplexMatrix(_N)
    print("\n\nrandom matrix =\n",np.asmatrix(random2))
    radius_array = gershGorinRadiusArray(random2)
    print("\n\napplying gershgorins theorem.")
    print("\n\nThe radii are ",radius_array)
    for i in range(_N):
        print("Eigenvalue ",i+1," : ",random2[i][i]," , radius = ",radius_array[i])


    # --- e) ---- plot convergence of matrix in semilogarithmic plot from jacobi algorithm
    print("\n\n--- e) ---------------------")
    plt.figure()
    plt.semilogy(np.arange(len(S)),S)           # plot semilogarithmic
    plt.title("")
    plt.grid()
    plt.xlabel("iterations of rotation")
    plt.ylabel("Convergence S")
    plt.title("Convergence of matrix along rotation")
    plt.show()


    # --- f) --- repeated jacobi transformations with N = 3 .. 50, plot time as function of N--------------------------
    print("\n\n--- f) ---------------------")
    Ns = np.arange(3,31)
    times = np.zeros(len(Ns), dtype= float)
    print("starting repeating jacobi transformation with NxN = 3x3 ... 50x50")
    time1 = time.time()
    for i in range(len(Ns)):                                # do rotations of all matrices
        A = randomMatrix(Ns[i])                             # random matrix
        print("N is ",str(i+3))
        [D, U, count, S,rads] = jacobi(A, 1e-1,1e5,0.5)     # jacobi method
        times[i] = time.time()-time1                       # measure time
        time1 = time.time()

    # plot size of matrix against rotation time
    plt.figure()
    plt.plot(Ns,times)
    plt.title("Process time t(N)")
    plt.grid()
    plt.xlabel("N x N Matrix")
    plt.ylabel("time t")
    plt.show()

    return

if __name__ == '__main__':
    main()