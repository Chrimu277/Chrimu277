
# Normal modes: vibrations, eigenfunctions of Translation operator
    #u_x lmn = u_x k *exp ( i(lka+mka+nka-wt) )
    # fcc or bcc --> make equation m*d2/dt2 u_x lmn = ..... (netwonians equation of motion with linear springs)
    # then converts equations into algebraics .. matrix !
    #

import numpy as np
from enum import Enum
import latticegeometry as LG
import matplotlib.pyplot as plt
from sympy import *

def u_lmn_bcc_vec(uk:np.ndarray,l:int,m:int,n:int,k:np.ndarray,a1:np.ndarray,a2:np.ndarray,a3:np.ndarray):
    return uk * np.exp(complex(0,1)*(l*np.dot(k,a1)+m*np.dot(k,a2)+n*np.dot(k,a3)))

def u_lmn_bcc_scal(uk:np.ndarray,l:int,m:int,n:int,kx:float,ky:float,kz:float,a:float):
    return uk * np.exp(complex(0,1)*(l+m)*kx*a/2)*np.exp(complex(0,1)*(m+n)*ky*a/2)*np.exp(complex(0,1)*(l+n)*kz*a/2)

# from some student project
# TODO WRONG ?
def algebraic_BCC_matrix(kx, ky, kz, a):
    # line 1
    m11 =  4 - np.cos(0.5*a*(kx+ky+kz)) - np.cos(0.5*a*(3*kx-ky-kz))\
           - np.cos(0.5*a*(-kx+3*ky-kz)) - np.cos(0.5*a*(-kx-ky+3*kz))

    m12 = - np.cos(0.5*a*(kx+ky+kz)) - np.cos(0.5*a*(3*kx-ky-kz))\
           + np.cos(0.5*a*(-kx+3*ky-kz)) - np.cos(0.5*a*(-kx-ky+3*kz))

    m13 = - np.cos(0.5*a*(kx+ky+kz)) + np.cos(0.5*a*(3*kx-ky-kz))\
           - np.cos(0.5*a*(-kx+3*ky-kz)) + np.cos(0.5*a*(-kx-ky+3*kz))
    # line 2
    m21 = m12
    m22 = m11
    m23 = - np.cos(0.5*a*(kx+ky+kz)) - np.cos(0.5*a*(3*kx-ky-kz))\
           + np.cos(0.5*a*(-kx+3*ky-kz)) + np.cos(0.5*a*(-kx-ky+3*kz))

    m31 = m13
    m32 = m23
    m33 = m11

    M = np.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
    return M

# http://lampx.tugraz.at/~hadley/ss1/phonons/bcc/bcc.php
def BCC_matrix_hadley1(kx:float, ky:float, kz:float, a:float, C21:float = 0):
    m11 = (2 / 3) * (-np.cos(0.5 * (-kx*a + ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a - kz*a)) - np.cos(
        0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)) + 4) - 2 * C21 * (np.cos(kx*a) - 1)
    m12 = (2 / 3) * (np.cos(0.5 * (-kx*a + ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a - kz*a)) + np.cos(
        0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)))
    m13 = (2 / 3) * (np.cos(0.5 * (-kx*a + ky*a + kz*a)) + np.cos(0.5 * (kx*a + ky*a - kz*a)) - np.cos(
        0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)))
    #m22 = (2 / 3) * (-np.cos(0.5 * (-kx*a + ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a - kz*a)) - np.cos(
    #    0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)) + 4) - 2 * C21 * (np.cos(ky*a) - 1)
    m23 = (2 / 3) * (-np.cos(0.5 * (-kx*a + ky*a + kz*a)) + np.cos(0.5 * (kx*a + ky*a - kz*a)) + np.cos(
        0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)))
    #m33 = (2 / 3) * (-np.cos(0.5 * (-kx*a + ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a - kz*a)) - np.cos(
    #    0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)) + 4) - 2 * C21 * (np.cos(kz*a) - 1)

    m21 = m12
    m22 = (2 / 3) * (-np.cos(0.5 * (-kx*a + ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a - kz*a)) - np.cos(
        0.5 * (kx*a - ky*a + kz*a)) - np.cos(0.5 * (kx*a + ky*a + kz*a)) + 4) - 2 * C21 * (np.cos(ky*a) - 1)

    m31 = m13
    m32 = m23
    m33 = (2 / 3) * (-np.cos(0.5 * (-kx * a + ky * a + kz * a)) - np.cos(0.5 * (kx * a + ky * a - kz * a)) - np.cos(
        0.5 * (kx * a - ky * a + kz * a)) - np.cos(0.5 * (kx * a + ky * a + kz * a)) + 4) - 2 * C21 * (
                      np.cos(kz * a) - 1)

    M = np.matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    return M

# http://lampx.tugraz.at/~hadley/ss1/phonons/fcc/fcc.php
# http://lampx.tugraz.at/~hadley/ss1/phonons/fcc/fcc2.php
def FCC_matrix_hadley1(kx:float, ky:float, kz:float, a:float, C21:float = 0):
    m11 = 4 - np.cos(0.5*a*(kx+ky)) - np.cos(0.5*a*(kx+kz)) - np.cos(0.5*a*(kx-ky)) - np.cos(0.5*a*(kx-kz))
    m22 = 4 - np.cos(0.5*a*(kx+ky)) - np.cos(0.5*a*(ky+kz)) - np.cos(0.5*a*(kx-ky)) - np.cos(0.5*a*(ky-kz))
    m33 = 4 - np.cos(0.5*a*(kx+kz)) - np.cos(0.5*a*(ky+kz)) - np.cos(0.5*a*(kx-kz)) - np.cos(0.5*a*(ky-kz))

    m12 = - np.cos(0.5*a*(kx+ky)) + np.cos(0.5*a*(kx-ky))
    m13 = - np.cos(0.5*a*(kx+kz)) + np.cos(0.5*a*(kx-kz))
    m23 = - np.cos(0.5*a*(ky+kz)) + np.cos(0.5*a*(ky-kz))

    m31 = m13
    m32 = m23
    m21 = m12

    M = np.matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    return M


def bcc_lambda(m,w,c):
    lambd = m*w**2 /(np.sqrt(3)*c)
    return np.matrix([lambd,0,0],[0,lambd,0],[0,0,lambd])

# Class to retrieve U V W factors for b1,b2,b3 for BCC symmetry points
# k = u b1 + v b2 + w b3
class SymmetryPoints_BCC(object):
    GAMMA = [0, 0, 0]
    H = [-0.5, 0.5, 0.5]
    P = [0.25, 0.25, 0.25]
    N = [0, 0.5, 0]

    def __init__(self):
        pass

    def getUVW(self, name:str) -> list:
        ret = None
        if name == "G" or name == "g" or name == "gamma" or name == "GAMMA":
            ret = self.GAMMA
        elif name == "H" or name == "h":
            ret = self.H
        elif name == "p" or name == "P":
            ret = self.P
        elif name == "n" or name == "N":
            ret = self.N
        else:
            print("dont know that symmetry point")
            ret = [-1,-1,-1]
        return ret

class SymmetryPoints_FCC(object):
    GAMMA = [0, 0, 0]
    X = [0, 1/2, 1/2]
    L = [1/2, 1/2, 1/2]
    W = [1/4, 3/4, 1/2]
    U =[1/4, 5/8, 5/8]
    K = [1/4, 3/4, 3/8]

    def __init__(self):
        pass

    def getUVW(self, name:str) -> list:
        ret = None
        if name == "G" or name == "g" or name == "gamma" or name == "GAMMA":
            ret = self.GAMMA
        elif name == "X" or name == "x":
            ret = self.X
        elif name == "L" or name == "l":
            ret = self.L
        elif name == "W" or name == "w":
            ret = self.W
        elif name == "U" or name == "u":
            ret = self.U
        elif name == "K" or name == "k":
            ret = self.K
        else:
            print("dont know that symmetry point")
            ret = [-1,-1,-1]
        return ret

def testBCCdispersion(m, a, c1=1, c21=0.1, N_interp=100, points = ["g", "h", "p", "g", "n"] ):
    # get dispersion relation of iron for Gamma to H
    a1_bcc, a2_bcc, a3_bcc = LG.primitiveLatticeVectorsBCC(a)  # primitive lattice vectors
    b1_bcc, b2_bcc, b3_bcc = LG.reciprocalVectors(a1_bcc, a2_bcc, a3_bcc)  # reciprocal lattice vectors

    # get symmetry points
    SPBCC = SymmetryPoints_BCC()                              # select symmetry points
    uvws = [SPBCC.getUVW(points[i]) for i in range(len(points))]    # get uvw coordinates for all SP's
    kxyzs = [(uvws[i][0] * b1_bcc + uvws[i][1] * b2_bcc + uvws[i][2] * b3_bcc) # get k-vector coordinates for all SP's
             for i in range(len(uvws))]

    # make k-vector arrays for distance between Symmetry points
    k_all, N_interp_is = ks_and_interpolation(kxyzs, N_interp)

    # solve
    energies, ks = solveAlgebraicEquations(k_all, N_interp_is, a, c21, BCC_matrix_hadley1)

    # plot EV over index of k-state
    plt.plot(np.arange(len(energies[0])), energies[0])
    plt.plot(np.arange(len(energies[1])), energies[1])
    plt.plot(np.arange(len(energies[2])), energies[2])
    maxenergy = np.matrix.max(np.asmatrix(energies))
    for j in range(len(points)):
        plt.vlines(x=np.sum(N_interp_is[:j]), ymin=0, ymax=10, colors="red", label=points[j])
        plt.text(x=np.sum(N_interp_is[:j]), y=-0.3, s=points[j])

    plt.title(f"Phonon Dispersion Relation of BCC Structure\n"
              f"{points}; C21 = {c21}")
    plt.grid()
    plt.ylabel("w * sqrt(m/c)")
    plt.xlabel("index of k-vector")
    plt.ylim((-0.5, 1.5 * maxenergy))
    plt.savefig("bcc_dispersion.jpg")
    plt.show()

def testFCCdispersion(m, a, c1=1, c21=0.1, N_interp=100, points = ["g", "x", "w", "k", "g","l"]):
    # get dispersion relation of iron for Gamma to H
    a1_fcc, a2_fcc, a3_fcc = LG.primitiveLatticeVectorsFCC(a)  # primitive lattice vectors
    b1_fcc, b2_fcc, b3_fcc = LG.reciprocalVectors(a1_fcc, a2_fcc, a3_fcc)  # reciprocal lattice vectors

    # get symmetry points
    SPFCC = SymmetryPoints_FCC()
    uvws = [SPFCC.getUVW(points[i]) for i in range(len(points))]  # get uvw coordinates for all SP's
    kxyzs = [(uvws[i][0] * b1_fcc + uvws[i][1] * b2_fcc + uvws[i][2] * b3_fcc)  # get k-vector coordinates for all SP's
             for i in range(len(uvws))]

    # make k-vector arrays for distance between Symmetry points
    k_all,N_interp_is = ks_and_interpolation(kxyzs, N_interp)

    # solve problem
    energies, ks = solveAlgebraicEquations(k_all, N_interp_is, a, c21, FCC_matrix_hadley1)

    maxenergy = np.max(energies)
    plt.figure()
    plt.plot(np.arange(len(energies[0])), energies[0])
    plt.plot(np.arange(len(energies[1])), energies[1])
    plt.plot(np.arange(len(energies[2])), energies[2])
    for j in range(len(points)):
        plt.vlines(x=np.sum(N_interp_is[:j]), ymin=0, ymax=10, colors="red", label=points[j])
        plt.text(x=np.sum(N_interp_is[:j]), y=-0.3, s=points[j])

    plt.title(f"Phonon Dispersion Relation of FCC Structure\n"
              f"{points}; C21 = {c21}")
    plt.grid()
    plt.ylabel("w * sqrt(m/c)")
    plt.xlabel("index of k-vector")
    plt.ylim((-0.5, 1.5 * maxenergy))
    plt.savefig("fcc_dispersion2.jpg")
    plt.show()

def calculatePhononDispersion(m, a, c1=1, c21=0.1, N_interp=100,
    matrixFunc=FCC_matrix_hadley1, SPclass=SymmetryPoints_FCC, points = ["g", "x", "w", "k", "g","l"]):
    # get dispersion relation of iron for Gamma to H
    a1_fcc, a2_fcc, a3_fcc = LG.primitiveLatticeVectorsFCC(a)  # primitive lattice vectors
    b1_fcc, b2_fcc, b3_fcc = LG.reciprocalVectors(a1_fcc, a2_fcc, a3_fcc)  # reciprocal lattice vectors

    # get symmetry points
    SP = SPclass()
    uvws = [SP.getUVW(points[i]) for i in range(len(points))]  # get uvw coordinates for all SP's
    kxyzs = [(uvws[i][0] * b1_fcc + uvws[i][1] * b2_fcc + uvws[i][2] * b3_fcc)  # get k-vector coordinates for all SP's
             for i in range(len(uvws))]

    # make k-vector arrays for distance between Symmetry points
    k_all,N_interp_is = ks_and_interpolation(kxyzs, N_interp)

    # solve problem
    energies, ks = solveAlgebraicEquations(k_all, N_interp_is, a, c21, matrixFunc)

    return energies, ks, N_interp_is

def plotPhononDispersion(energies, ks, N_ints, points, title):
    # plot EV over index of k-state
    plt.scatter(np.arange(len(energies[0])), energies[0])
    plt.scatter(np.arange(len(energies[1])), energies[1])
    plt.scatter(np.arange(len(energies[2])), energies[2])
    maxenergy = np.matrix.max(np.asmatrix(energies))
    for j in range(len(points)):
        plt.vlines(x=np.sum(N_ints[:j]), ymin=0, ymax=10, colors="red", label=points[j])
        plt.text(x=np.sum(N_ints[:j]), y=-0.3, s=points[j])

    plt.title(f"{title}\n"
              f"{points}")
    plt.grid()
    plt.ylabel("w * sqrt(m/c)")
    plt.xlabel("index of k-vector")
    plt.ylim((-0.5, 1.5 * maxenergy))
    plt.show()
# "BCC" or "FCC"
def solveAlgebraicEquations(k_all, N_interp_is, a, c21, matrixFunction = FCC_matrix_hadley1):

    lambdas = []  # stores eigenvalues
    ks = []  # stores |k|'s
    plt.figure()
    for j in range(len(k_all)):  # for each SymPoint <--> SymPoint distance
        for i in range(N_interp_is[j]):  # for each interpolation step (k-vector)
            # get eigenvalue .. lambda = m w² / sqrt3 c
            # make equation M u_xyz = Mw² / sqrt3 C  * u_xyz
            # M = algebraic_BCC_matrix(k_all[j][0][i], k_all[j][1][i], k_all[j][2][i], a=a_bccIron)
            M = matrixFunction(k_all[j][0][i], k_all[j][1][i], k_all[j][2][i], a=a, C21=c21)

            # solve eigenvalue equation M(v-IL) = 0
            eval, evec = np.linalg.eig(M)
            # go through eigenvalues
            eval = [complex(eval[i]).real for i in range(len(eval))]
            for i in range(len(eval)):
                if eval[i] < 0:
                    eval[i] = 0
                eval[i] = np.sqrt(eval[i].real)  # *np.sqrt(3)
                # append eigenvalues and k-vectors
            lambdas.append(np.sort(np.array(eval)))
                # append absolute values of k
            ks.append(k_all[j][0][i] ** 2 + k_all[j][1][i] ** 2 + k_all[j][2][i] ** 2)

    # transpose the resulting eigenvalues for plotting
    lambdas2 = np.asarray(lambdas).transpose()
    return lambdas2,ks

# takes vector of [kx,ky,kz]_i and interpolates them
def ks_and_interpolation(kxyzs,N_interp):
    k_all = []
    d0 = kxyzs[1] - kxyzs[0]
    abs_d0 = np.sqrt(np.dot(d0, d0))
    N_interp_is = []
    for i in range(len(kxyzs) - 1):
        di = kxyzs[i + 1] - kxyzs[i]
        abs_di = np.sqrt(np.dot(di, di))
        N_interp_i = int(N_interp * abs_di / abs_d0)
        k_all.append([np.linspace(kxyzs[i][0], kxyzs[i + 1][0], N_interp_i),
                      np.linspace(kxyzs[i][1], kxyzs[i + 1][1], N_interp_i),
                      np.linspace(kxyzs[i][2], kxyzs[i + 1][2], N_interp_i)])
        N_interp_is.append(N_interp_i)

    return k_all,N_interp_is

if __name__ == '__main__':

    # SET SIMULATION PARAMETERS
    amu = 1.66054e-27       # atomic mass unit
    m_Carbon = 12*amu       # atomic mass of Carbon C12
    m_Iron = 55.845*amu     # atomic mass of Iron , Fe
    a_bccIron = 0.2866e-9   # meters, Lattice constant for BCC Iron
    a_fccIron = 0.3571e-9   # Lattice constant for FCC Iron

    w = 1e+12               # vibration frequency (NOT NEEDED)
    c1 = 1e0                 # spring constant nearest neighbours
    c21 = 0.3*c1            # spring constant next nearest neighbours

    points_bcc = ["g", "h", "p", "g", "n"]
    points_fcc = ["g", "x", "w", "k", "g","l"]
    N_interpolate = 100

    #e_bcc, k_bcc, N_int_bcc = calculatePhononDispersion(m=m_Iron,a=a_bccIron,N_interp=N_interpolate, c1=c1,c21=c21,
    #matrixFunc = BCC_matrix_hadley1,SPclass = SymmetryPoints_BCC, points = ["g", "h", "p", "g", "n"])

    #plotPhononDispersion(e_bcc,k_bcc,N_int_bcc,points_bcc,"BCC thing")
    #e_fcc, k_fcc, N_int_fcc = calculatePhononDispersion(m=m_Iron,a=a_fccIron, N_interp=N_interpolate, c1=c1, c21 = c21,
    #matrixFunc = FCC_matrix_hadley1, SPclass = SymmetryPoints_FCC, points = ["g", "x", "w", "k", "g","l"])

    testBCCdispersion(m=m_Carbon,a=a_bccIron,c1=c1,c21=c21,N_interp=N_interpolate,points=points_bcc)
    testFCCdispersion(m_Carbon,a_fccIron,c1,c21,N_interpolate,points_fcc)


