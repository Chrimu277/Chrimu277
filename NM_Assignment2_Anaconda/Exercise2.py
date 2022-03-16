# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt


# ----FUNCTIONS DEFINITIONS--------
def y_fda1(xi,h):
    return 0

# analytical solution from assignment sheet
def y_analytical(x,d,k,y0,yd): # 1/100 im argument, da sons werte zu groß
    y = (yd*np.sinh(k*x) + y0*np.sinh(k*(d-x))) / np.sinh(k*d)
    return y;

# derivate of analytical solution y(x)
def y_analytical_derivative1(x,d,k,y0,yd): # 1/100 im argument, da sons werte zu groß
    y_strich = k*(yd*np.cosh(k*x) - y0*np.cosh(k*(d-x))) / np.sinh(k*d)
    return y_strich;

# uses the analytical derivative to calculate the charge on the capacitor
def charge(y0,yd,d):
    y_strich_0 = y_analytical_derivative1(0,d,_K,y0,yd)
    y_strich_d = y_analytical_derivative1(d,d,_K,y0,yd)
    return (y_strich_0 - y_strich_d) * 2*_C0/(_K**2)

# uses the formula for numerical derivatives on the numerical solution to calculate charge
def charge_finit_element_numerical(y0,yd,d, n = 1000):
    y_xi = get_yxi_numerical(n,d,y0,yd,_K);
    y_strich_0 = (y_xi[1]-y_xi[0]) / (d/n)
    y_strich_d = (y_xi[n-1] - y_xi[n-2]) / (d/n)

    return (y_strich_0 - y_strich_d) * 2*_C0/(_K**2)

# uses the formula for numerical derivatives on the analytical solution to calculate charge
def charge_finit_element_analytical(y0,yd,d, fac = 1000):

    y_strich_0 = (y_analytical((d/fac),d,_K,y0,yd) - y_analytical(0,d,_K,y0,yd)) / (d/fac)
    y_strich_d = (y_analytical(d,d,_K,y0,yd) - y_analytical(d-(d/fac),d,_K,y0,yd)) / (d/fac)

    return (y_strich_0 - y_strich_d) * 2*_C0/(_K**2)


# n: length of return vector yxi[n]
# d: capacitor width
# y0,yd: randbedignungen
# k: debye length
# calculates y(xi) by inverse matrix multiplication on b; x = A_in * b
def get_yxi_numerical(n,d,y0,yd,k):
    # create result-vector b
    b = np.zeros(n)
    b[0] = y0
    b[n - 1] = yd

    # create and fill matrix A
    A = [[0 for i in range(n)] for j in range(n)]  # A is m x n .. m rows, n cols

    for i in range(n):
        for j in range(n):
            if i == j:
                A[i][j] = 2 + (d / n * k) ** 2
            elif i == j + 1:
                A[i][j] = -1
            elif j == i + 1:
                A[i][j] = -1

    # invert A to calculate y(xi)
    A_in = np.linalg.inv(A)

    #y = A^-1 * b
    y_xi_found = A_in.dot(b)

    return y_xi_found


# --------VARIABLES -----------

# boundary conditions
_Y0 = -0.5
_YD = 0.5

# constants
_K = 0.323738

_H = 1           # step , xi = h*i
_C0 = 0.006     # 1/nm^3

_D = 100        # nm
_N = 1000       # steps for 0...d

# ------START -------------------

# ----c) solve 1d debye hückel equation

# x element [0,d]
_Xis = np.linspace(0, _D, _N)      # nm   xi [0..100]

print("N = "+str(_N))
print("len Xis = "+str(len(_Xis)))

# calculate array of solution both numerical , analytical
y_xi_found = get_yxi_numerical(_N,_D,_Y0,_YD,_K)
y_xi_analytical = y_analytical(_Xis, _D, _K, _Y0, _YD)

print("len yxi found "+str(len(y_xi_found)))
print("len yxi analytical "+str(len(y_xi_analytical)))

# plot of both y(xi)
plt.figure()
plt.plot(_Xis, y_xi_analytical, 'r', label = 'analytical')
plt.title("Analytical y(xi)")
plt.grid()
plt.show()
plt.plot(_Xis, y_xi_found, 'b', label = 'found')
plt.title("Matrix-calculated y(xi)")
plt.grid()
plt.show()

# --------d) implement script to calculate capacitance of supercapaicator
# C = dq / dy0  = q1 - q2  /  y0_1  - y0_2
# q = 2 c0 / k^2  * (y'(0) - y'(d)

# randbedingungen
Y01 = -0.3
YD1 = 0.5

Y02 = -0.4
YD2 = 0.5

# create array of d's
Ds = np.linspace(1,100,100) # dy = 1

#---method 1 to calc charge
# qs11 :float = []
# # qs21 :float= []
# # for i in range(len(Ds)):
# #     qs11.append(charge_finit_element_numerical(Y01,YD1,Ds[i]))
# #     qs21.append(charge_finit_element_numerical(Y02,YD2,Ds[i]))
# #
# # Cs1 = []
# # for i in range(len(Ds)):
# #     Cs1.append((qs11[i]-qs21[i])/(Y01-Y02))
# #
# # plt.xlabel("d / nm")
# # plt.ylabel("q")
# # plt.title("Charge q method 1")
# # plt.plot(Ds,qs11,'b',label = 'Q1i(d)')
# # plt.plot(Ds,qs21,'g',label = 'Q2i(d)')
# # plt.show()
# #
# # plt.xlabel("d / nm")
# # plt.ylabel("Ci(d) / F ")
# # plt.plot(Ds,Cs1,'r')
# # plt.title("Capacitance with method 1")
# # plt.show()

# ---method 2 to calc charge
qs12 :float = charge(Y01,YD1,Ds)
qs22 :float = charge(Y02,YD2,Ds)
Cs2 = (qs22-qs12) / (Y02-Y01)

plt.xlabel("d / nm")
plt.ylabel("q")
plt.title("Charge q method 2")
plt.plot(Ds,qs12,'b',label = 'Q1(d)')
plt.plot(Ds,qs22,'g',label = 'Q2(d)')
plt.show()
plt.xlabel("d / nm")
plt.ylabel("C(d) / F ")
plt.plot(Ds,Cs2,'r')
plt.title("Capacitance with method 2")
plt.show()

# ---method 3 to calc charge

qs13 = charge_finit_element_analytical(Y01,YD1,Ds)
qs23 = charge_finit_element_analytical(Y02,YD2,Ds)
Cs3 = (qs13-qs23) / (Y01-Y02)

plt.xlabel("d / nm")
plt.ylabel("q")
plt.title("Charge q method 3")
plt.plot(Ds,qs13,'b',label = 'Q1(d)')
plt.plot(Ds,qs23,'g',label = 'Q2(d)')
plt.show()
plt.xlabel("d / nm")
plt.ylabel("C(d) / F ")
plt.plot(Ds,Cs3,'r')
plt.title("Capacitance with analytical derivative")
plt.show()