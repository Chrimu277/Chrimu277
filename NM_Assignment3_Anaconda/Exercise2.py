# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.linalg import fractional_matrix_power


# decompose pixel matrix into 3 rgb matrices
def RGB_Image(pixels):
    r = np.zeros(shape=(len(pixels),len(pixels[0])))
    g = r.copy()
    b = r.copy()
    for row in range(len(pixels)):
        for col in range(len(pixels[0])):
            r[row][col] = pixels[row][col][0]
            g[row][col] = pixels[row][col][1]
            b[row][col] = pixels[row][col][2]
    return r,g,b


# checks properties of matrix
def checkMatrix(matrix):

    print("dimensions: ",np.shape(matrix))

    # check symmetry
    symmetricBool = checkSymmetry(matrix)
    print("matrix symmetric: ", symmetricBool)

    # check unitary
    unitaryBool = checkUnitary(matrix)
    print("matrix unitary: ", unitaryBool)

    # check diagonalizability
    diagonalizableBool = checkDiagonalizable(matrix)
    print("matrix diag.-able: ", diagonalizableBool)

    # check diagonal dominance
    diagDomBool = checkDiagDominance(matrix)
    print("matrix diag dom..: ", diagDomBool)

    # check sparsity
    sparseBool = checkSparse(matrix)
    print("matrix sparse: ", sparseBool)

    # check regularity
    regularBool = checkRegular(matrix)
    print("matrix regular: ", regularBool)
    return


# takes only a square matrix of topleft corner
def upperSquare(pixels):
    shape = np.shape(pixels)
    sidelength = shape[0] if shape[0] < shape[1] else shape[1]
    usqm = [[0 for i in range(sidelength)] for j in range(sidelength)]

    for row in range(sidelength):
        for col in range(sidelength):
            usqm[row][col] = pixels[row][col]

    return usqm


# checks matrix for symmetry
def checkSymmetry(matrix, tol = 1e-9):
    # u - ut = 0
    matrix_2 = matrix - np.transpose(matrix)
    return checkMatrixZero(matrix_2)


# checks is a matrix is unitary (Ut = U^-1)
def checkUnitary(matrix, tol = 1e-9):
    # Ut = u^-1
    inv = np.linalg.inv(matrix)
    shouldZero = inv-np.transpose(matrix).conjugate()
    return checkMatrixZero(shouldZero)


#checks if a matrix is (almost) 0 in all elements
def checkMatrixZero(matrix, tol = 1e-9):
    AbsSum = 0.0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            AbsSum += np.abs(matrix[i][j])
    return True if AbsSum < tol else False


#checks if a matrix is "normal" in other words "diagonalizable"
def checkDiagonalizable(matrix):
    AAt = np.dot(matrix, np.transpose(matrix).conjugate())
    AtA = np.dot(np.transpose(matrix).conjugate(), matrix)
    return checkMatrixZero((AAt-AtA))   # true if AAt = AtA


# get absolute sum of all off diag. elements.
def getSumNebenDiagElemente(A):
    A_temp = A.copy();
    S = 0.0
    for i in range(len(A_temp)):
        for j in range(len(A_temp)):
            if(i != j):
                S += np.abs(A_temp[i][j])
    return S;


# dominant diag if every diagonal element is biggest than ALL others
def checkDiagDominance(matrix):
    matrix_temp = matrix.copy()
    sum_others = getSumNebenDiagElemente(matrix_temp) # gets the sum of the absolutes of the offdiag elements
    diagDom = True

    for i in range(len(matrix_temp)):
        if matrix_temp[i][i] < sum_others:  # if every diagonal element is bigger than the sum of the offdiagonal elements,
            diagDom = False;                # its diagonally dominant

    return diagDom


# checks matrix for sparsity
def checkSparse(matrix, tol = 1e-9):
    matrix_temp = matrix.copy()
    n = len(matrix_temp)
    nr_elements = n**2
    nr_zeroes = 0;

    #count the empty places of the matrix
    for i in range(n):
        for j in range(n):
            if matrix_temp[i][j] < tol:     # if morse than 50% of matrix is empty, its sparse
                nr_zeroes += 1;

    sparsity = nr_zeroes/nr_elements        # factor between 0 ... 1

    return True if sparsity > 0.5 else False


# check if matrix is regular
def checkRegular(matrix):
    matrix_temp = matrix.copy()
    regular = True
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix_temp[i][j] <= 0:  # if any element is equal or smaller than 0, false
                regular = False

    return regular;


# powermethod to estimate the biggest eigenvalue and corresponding vektor
def powerMethod(matrix, maxSteps = 1e5, tol = 1e-5):

    A = matrix.copy()   # copy out matrix

    RESIDUA = []            # residuen
    X_N = np.random.rand(len(A)) % 1 # eigenvektor
    LAMBDA = 0              # eigenwert


    #find largest eigenvektor, eigenvalue
    for i in range(int(maxSteps)):

        # x_p+1 = Ax_p
        x_p = np.dot(A,X_N)
        x_p = x_p / (np.linalg.norm(x_p))    # norm

        # EW = Xn_T * M * Xn / (Xn_T Xn)
        lambda_p = np.dot(x_p.transpose().conjugate(), np.dot(A, x_p)) / (np.dot(x_p.transpose().conjugate(), x_p))

        # residuum
        RESIDUA.append(lambda_p)

        X_N = x_p  # set eigenvector
        if np.abs(LAMBDA-lambda_p) < tol:
            break;
        LAMBDA = lambda_p   # set eigenwert

    for i in range(len(RESIDUA)):
        RESIDUA[i] = np.abs(RESIDUA[i] - LAMBDA)

    #now x_n is largest eigenvalue
    return X_N, LAMBDA, RESIDUA


# gets the n largest eigenvalues and eigenvectors of B
def deflation(matrix, n):
    B = matrix.copy()

    EW = []    # eigenvektors, eigenvalues, residua
    EV = []
    RES = []

    for i in range(n):

        ev, ew, RES = powerMethod(B,1e3,1e-5)         # power method
        B = B - ew* np.outer(ev,ev.transpose().conjugate())  # deflation of matrix
        EW.append(ew)
        EV.append(ev)
        RES.append(RES)

    return EV,EW,RES

# gram schmidt orthogonalization to calculate a orthogonalized set of base vectors u[][]
def gramSchmidt(u):

    n = len(u)  # count of vectors
    i = len(u[0])   # length of 1 vector

    base = [[0 for i in range(i)] for j in range(n)]

    # get all vk
    for k in range(n):
        s = 0.0
        for j in range(k):
                s +=  base[j] * ((u[k].dot(base[j])) / (base[j].dot(base[j])))

        vk = u[k] - s                       # gram schmid formula
        vk = vk / (np.linalg.norm(vk))      # norm the vektor
        base[k] = vk.copy()
    return base


# generates an SVD reduced image of 3 color matrices (r,g,b) and saves it under filename
def makeImage(red,green,blue,n):
    filename = "Vitruvian"+str(n)+".bmp"

    # make SVD representation of red, green, blue
    redSVD = np.transpose(singValDecomp(red,n))
    greenSVD = np.transpose(singValDecomp(green, n))
    blueSVD = np.transpose(singValDecomp(blue, n))

    #redSVD = np.zeros(np.shape(redSVD))
    #greenSVD = np.zeros(np.shape(redSVD))
    #blueSVD = np.zeros(np.shape(redSVD))

    # set together pixels
    pixels = [[ [redSVD[i][j], greenSVD[i][j], blueSVD[i][j]]  for i in range(len(redSVD))] for j in range(len(redSVD[0]))]
    image = Image.fromarray(np.uint8(pixels))
    # save image
    image.save(filename)
    return

# function for the n-SingularValueDecomposition of a matrix (n repeated deflatio)
def singValDecomp(matrix,n):
    #assignment: SVD
    A = matrix.copy()
    B = np.dot(np.transpose(A).conjugate(),A)
    evs,ews,residua = deflation(B, n)         # get n biggest eigenvalues and eigenvectors
    V = gramSchmidt(evs)                        # V is orthonormalized egenvectors
    D = np.zeros((len(ews),len(ews)))

    for i in range(len(ews)):
        D[i][i] = ews[i]

    #L = D^(1/2)
    Lambda = fractional_matrix_power(D, 0.5)

    #Un = A V D^-1/2
    Un = np.dot(A, np.dot(np.transpose(V).conjugate(), (fractional_matrix_power(D, -0.5))))

    #An = Un*L*Vt
    An = np.dot(Un, (np.dot(Lambda, (np.conjugate(V)))))
    An = np.absolute(An)    # 0 .. x > 0
    max = np.max(An)        # norm An auf 255
    An = An*(max/255)
    return An

# ----------------------------------------------------------------------------------------------------------------------
# START
# ----------------------------------------------------------------------------------------------------------------------

def main():
    # open image
    print("\nopening image")
    im = Image.open("Vitruvian.bmp")
    print("\nim is type: ",type(im))

    # make im to pixelmatrix
    pixels = np.asarray(im)
    print("\ndimensions are ",np.shape(pixels))

    print("\ntaking top left corner:")
    # color decompose .. these are the  Av
    r_matrix, g_matrix, b_matrix = RGB_Image(pixels)
    print("\nits dimensions are ",np.shape(r_matrix))

    # show image full
    plt.figure()
    plt.title("Original File Data")
    plt.imshow(im)
    plt.show()

    # show image redvalues
    upSq_red = upperSquare(r_matrix)
    plt.figure()
    plt.title("Original File Data")
    plt.imshow(upSq_red)
    plt.show()


    # --- a) --- check matrix properties ---------------------------------------------------------
    print("\n\n--- a) --- \n check upper square red matrix")
    checkMatrix(upSq_red)

    # --- b) --- check B = Avt Av
    print("\n\n--- b) --- \n check B")
    B = np.dot(np.transpose(r_matrix).conjugate(), r_matrix)
    print("B shape is ",B.shape)
    checkMatrix(B)

    # --- c) --- Power method ---------------------------------------------------------------------
    print("\n\n--- c) --- \n Power method")

    ev, ew, residua = powerMethod(B)
    print("\n\nLargest Eigenvektor (powermethod:\n ",ev, " \nCorresponding Eigenvalue: ", ew)
    print("\nresidua are: ",residua)
    #plot residua
    plt.figure()
    plt.plot(np.arange(len(residua)),residua)
    plt.title("Residua over iterations of power method")
    plt.xlabel("Iterations")
    plt.ylabel("Residua")
    plt.show()

    # --- d) --- Deflation trhough repeated von Mises method (power method)
    print("\n\n--- d) --- \n Deflation")
    n = 5
    eigenvalues, eigenvectors, residuas = deflation(B,n)
    print("\n\nDeflation: \nEigenvalues",eigenvalues)
    for i in range(len(eigenvectors)):
        print("ev",i," = ",eigenvectors[i]," .....")

    # --- e) --- perform a SVD of length n of all 3 colors to compress image   1
    print("\n starting n-singular value decomposition, with n= 1, 5, 10, 20")
    # for each n, make a nSVD composed image
    n_array = np.array([1,5,10,20])
    for i in range(len(n_array)):
        print("making image with n=", n_array[i])
        makeImage(r_matrix,g_matrix,b_matrix,n_array[i])


    print("finished")

if __name__ == '__main__':
    main();
