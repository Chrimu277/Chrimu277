# Muehleder Christoph 01604413
import numpy as np
import matplotlib.pyplot as plt

#Gauß Seidel Method for exercise b)
#now, V1,V2,V3,V4 = const, I1..4(t) are solved
def GS_b(steps):
    # return vector
    i1, i2, i3, i4 = [], [], [], []
    # matrix A
    a_ij = [[_R1,_R2,0,0],
            [0,-_R2,_R3,0],
            [0,0,_R3,0],
            [-1,1,1,1]]
    # A x = b
    b = [_V1-_V2,_V2-_V3,_V4-_V3,0]
    # starting vectors x
    xi_k = [0,0,0,0]
    xi_k_1 = [0,0,0,0]

    #formula for Gaußseidel method
    for c in range(steps):
        for i in range(4):
            sum1 = 0;
            sum2 = 0;
            for j in range(i):
                sum1+= a_ij[i][j]*xi_k_1[j]
            for j in range(i+1,4):
                sum2+= a_ij[i][j]*xi_k[j]
            xi_k_1[i] = (b[i] - sum1 - sum2) / a_ij[i][i]
        # old = new
        xi_k = xi_k_1
        # append to solution
        i1.append(xi_k_1[0])
        i2.append(xi_k_1[1])
        i3.append(xi_k_1[2])
        i4.append(xi_k_1[3])

    return i1,i2,i3,i4

#Gauß Seidel Method for exercise c)
#now, V1,V2,V3 = const, V4(t) and I1..4(t) are solved
def GS_c(time):
    # return vectors
    v1, v2, v3, v4, i1, i2, i3, i4 = [], [], [], [], [], [], [], []
    # Matrix A
    a_ij = [[_R1,_R2,0,0],
            [0,-_R2,_R3,0],
            [0,0,_R3,0],
            [-1,1,1,1]]
    # A x = b
    b = [_V1-_V2,_V2-_V3,0-_V3,0]
    volt4 = 0
    # starting vectors
    xi_k = [0,0,0,0]
    xi_k_1 = [0,0,0,0]

    #iterations fpr gauß seidel method
    for c in range(time):
        for i in range(4):
            sum1 = 0;
            sum2 = 0;
            for j in range(i):
                sum1+= a_ij[i][j]*xi_k_1[j]
            for j in range(i+1,4):
                sum2+= a_ij[i][j]*xi_k[j]
            xi_k_1[i] = (b[i] - sum1 - sum2) / a_ij[i][i]
        xi_k = xi_k_1

        i1.append(xi_k_1[0])
        i2.append(xi_k_1[1])
        i3.append(xi_k_1[2])
        i4.append(xi_k_1[3])
        v1.append(_V1)
        v2.append(_V2)
        v3.append(_V3)
        v4.append(volt4)

        volt4 += xi_k_1[3] / _C
        b[2] = volt4 - _V3

    return i1,i2,i3,i4,v1,v2,v3,v4

#Gauß Seidel Method for exercise d)
#now, V1(t) is given ( stimulus ) V2,V3 = const, V4(t) and I1..4(t) are solved
def GS_d(v1):
    time = len(v1)
    # return vectors
    v2, v3, v4, i1, i2, i3, i4 = [], [], [], [], [], [], []

    # matrix for calculation
    a_ij = [[_R1,_R2,0,0],
            [0,-_R2,_R3,0],
            [0,0,_R3,0],
            [-1,1,1,1]]
    # A x = b
    b = [_V1-_V2,_V2-_V3,_V4-_V3,0]
    # starting vectors
    xi_k = [0,0,0,0]
    xi_k_1 = [0,0,0,0]
    volt4 = 0;

    # gauß seidel formula
    for t in range(time):
        #calculate i-th component
        for i in range(4):
            sum1 = 0;
            sum2 = 0;
            for j in range(i):
                sum1+= a_ij[i][j]*xi_k_1[j]
            for j in range(i+1,4):
                sum2+= a_ij[i][j]*xi_k[j]
            xi_k_1[i] = (b[i] - sum1 - sum2) / a_ij[i][i]

        xi_k = xi_k_1

        #append to return vectors
        i1.append(xi_k_1[0])
        i2.append(xi_k_1[1])
        i3.append(xi_k_1[2])
        i4.append(xi_k_1[3])
        #v1 is argument, will be returned as is
        v2.append(_V2)
        v3.append(_V3)
        v4.append(volt4)

        #recalculate entries of b-vektor
        volt4 += xi_k[3] / _C
        b[2] = (volt4 - _V3)
        b[0] = (v1[t] - _V2) # new v1

    return i1,i2,i3,i4,v1,v2,v3,v4

# --------GLOBAL VARIABLES-----------#

# initialize variables
_V1 = 0
_V2 = 0.5
_V3 = -0.5
_V4 = 0

_R1 = 1
_R2 = 25
_R3 = 5
_C = 150

_i1, _i2, _i3, _i4, _v1, _v2, _v3, _v4 = [], [], [], [], [], [], [], []

# --------START-----------------------#

# -----c)----- Implement Gauss Seidel method for getting I_i, V_i
print("c)")

# calculate the resting currents
[_i1, _i2, _i3, _i4] = GS_b(5)

print("I1 = " + str(_i1[len(_i1) - 1]))
print("I2 = " + str(_i2[len(_i2) - 1]))
print("I3 = " + str(_i3[len(_i3) - 1]))
print("I4 = " + str(_i4[len(_i4) - 1]))

# ---- d)----- V4(t+dt) = V4(t) + I4(t) * dt/C
print("d)")

# calculate Vi(t) , Ii(t)

_Tmax = 1000
_dt = 1

_i1.clear()
_i2.clear()
_i3.clear()
_i4.clear()
_v1.clear()
_v2.clear()
_v3.clear()
_v4.clear()

#calculate all currents and v4(t)
[_i1, _i2, _i3, _i4, _v1, _v2, _v3, _v4] = GS_c(_Tmax);

#plot
plt.plot(np.arange(_Tmax), _i1, 'r', label ='I1')
plt.plot(np.arange(_Tmax), _i2, 'b', label ='I2')
plt.plot(np.arange(_Tmax), _i3, 'g', label ='I3')
plt.plot(np.arange(_Tmax), _i4, 'y', label ='I4')
plt.plot(np.arange(_Tmax), _v4, 'k', label ='V4')

plt.xlim(0, _Tmax)
plt.ylim(-0.5,0.5)

plt.title("Plot of d) Resting Potential")
plt.ylabel("I(t) / A   ...   V(t) / V")
plt.xlabel("t / s")

plt.grid('grey',linestyle='--', linewidth='1')
plt.legend()
plt.show();

print("V4(T) = " + str(_v4[len(_v4) - 1])) # resting potential


# ---- e) ---- stimulus action potential V1
print("e)")

_i1.clear()
_i2.clear()
_i3.clear()
_i4.clear()
_v1.clear()
_v2.clear()
_v3.clear()
_v4.clear()

#   make V1(t)
_V0 = 0.3
_TAU = 50
for t in range(_Tmax):
    _v1.append(_V0 * np.exp(-((t / _TAU) ** 2)))

#calculate all currents an voltages
[_i1, _i2, _i3, _i4, _v1, _v2, _v3, _v4] = GS_d(_v1)

# plot process
plt.plot(np.arange(_Tmax), _i1, 'r', label ='I1')
plt.plot(np.arange(_Tmax), _i2, 'b', label ='I2')
plt.plot(np.arange(_Tmax), _i3, 'g', label ='I3')
plt.plot(np.arange(_Tmax), _i4, 'y', label ='I4')
plt.plot(np.arange(_Tmax), _v1, 'k', label ='V1')
plt.plot(np.arange(_Tmax), _v4, label ='V4')

plt.xlim(0,_Tmax)
plt.ylim(-0.5,0.5)

plt.title("Plot of e) Action Potential")
plt.ylabel("I(t) / A   ...   V(t) / V")
plt.xlabel("t / s")
plt.grid('grey',linestyle='--', linewidth='1')

plt.legend()
plt.show()