# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt

#gauß seidel method, solve system of linera equations
def GS(steps):
    # return vector
    i1, i2, i3, i4, i5, i6 = [], [], [], [], [], []
    # matrix
    a_ij = [[_R1,   _R2,    0,      0,      0,      0],
            [1,     -1,     1,      0,      0,      0],
            [0,     0,      _R3,    _R4,    _R5,    0],
            [0,     0,      0,      _R4,    0,      _R6],
            [-1,     1,     0,      0,      -1,     0],
            [0,     0,      -1,     1,      0,      -1]]

    # a*x = b
    b = [_V1,0,_V2,_V3,0,0]

    # starting vectors
    xi_k = [0,0,0,0,0,0]
    xi_k_1 = [0,0,0,0,0,0]

    # formula for gauß seidel method
    for c in range(steps):
        for i in range(6):
            sum1 = 0;
            sum2 = 0;
            for j in range(i):
                sum1+= a_ij[i-1][j-1]*xi_k_1[j-1]
            for j in range(i+1,6):
                sum2+= a_ij[i-1][j-1]*xi_k[j-1]
            xi_k_1[i-1] = (b[i-1] - sum1 - sum2) / a_ij[i-1][i-1]

        # old vector = new vector
        xi_k = xi_k_1
        #append to solution
        i1.append(xi_k_1[0])
        i2.append(xi_k_1[1])
        i3.append(xi_k_1[2])
        i4.append(xi_k_1[3])
        i5.append(xi_k_1[4])
        i6.append(xi_k_1[5])

    return i1,i2,i3,i4,i5,i6

# succeessive over relaxation
def SOR(steps,omega):
    #return vectors
    i1, i2, i3, i4, i5, i6 = [], [], [], [], [], []
    #matrix
    a_ij = [[_R1,   _R2,    0,      0,      0,      0],
            [1,     -1,     1,      0,      0,      0],
            [0,     0,      _R3,    _R4,    _R5,    0],
            [0,     0,      0,      _R4,    0,      _R6],
            [-1,     1,     0,      0,      -1,     0],
            [0,     0,      -1,     1,      0,      -1]]

    #a * x = v
    b = [_V1,0,_V2,_V3,0,0]

    # starting vectors
    xi_k = [0,0,0,0,0,0]
    xi_k_1 = [0,0,0,0,0,0]

    #use formula for SOR method ( similar to GS)
    for c in range(steps):
        for i in range(6):
            sum1 = 0;
            sum2 = 0;
            for j in range(i):
                sum1+= a_ij[i][j]*xi_k_1[j]
            for j in range(i+1,6):
                sum2+= a_ij[i][j]*xi_k[j]

            xi_k_1[i] = (1-omega) * xi_k[i] + (b[i] - sum1 - sum2) * omega / a_ij[i][i]
        # old = new
        xi_k = xi_k_1
        # append to solution
        i1.append(xi_k_1[0])
        i2.append(xi_k_1[1])
        i3.append(xi_k_1[2])
        i4.append(xi_k_1[3])
        i5.append(xi_k_1[4])
        i6.append(xi_k_1[5])

    return i1,i2,i3,i4,i5,i6


# -----VARIABLE DEFINITIONS --------
_R1 = 2
_R2 = 4
_R3 = 1
_R4 = 2
_R5 = 2
_R6 = 4

_V1 = 10
_V2 = 17
_V3 = 14

i1, i2, i3, i4, i5, i6 = [], [], [], [], [], []

# --------START --------------

# ---calculate currents with GaußSeidel method---

_ITERATIONS_GS = 5

# get solution from GS method
[i1, i2, i3, i4, i5, i6] = GS(_ITERATIONS_GS);

# plot GS process
iterations = np.arange(_ITERATIONS_GS)
plt.plot(iterations,i1, label = 'i1')
plt.plot(iterations,i2, label = 'i2')
plt.plot(iterations,i3, label = 'i3')
plt.plot(iterations,i4, label = 'i4')
plt.plot(iterations,i5, label = 'i5')
plt.plot(iterations,i6, label = 'i6')
plt.legend()
plt.show()
#diverges

# ----calculate currents with successive over relaxation method---
_ITERATIONS_SOR = 10000;
_OMEGA_SOR = 0.01

# get solution with SOR method
[i1, i2, i3, i4, i5, i6] = SOR(_ITERATIONS_SOR,_OMEGA_SOR);

# plot sor process
iterations = np.arange(_ITERATIONS_SOR)
plt.plot(iterations,i1, label = 'i1')
plt.plot(iterations,i2, label = 'i2')
plt.plot(iterations,i3, label = 'i3')
plt.plot(iterations,i4, label = 'i4')
plt.plot(iterations,i5, label = 'i5')
plt.plot(iterations,i6, label = 'i6')
plt.legend()
plt.show()

print("i1...i6 = ")
print(i1[len(i1)-1])
print(i2[len(i1)-1])
print(i3[len(i1)-1])
print(i4[len(i1)-1])
print(i5[len(i1)-1])
print(i6[len(i1)-1])
