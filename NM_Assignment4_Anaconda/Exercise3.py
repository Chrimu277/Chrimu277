# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

# ----------------------------------------------------------------------------------------------------------------------

# the function for which X' = f(X,t,..)
def func1(x,t,m,q,B):

    f1 = x[3].copy()
    f2 = x[4].copy()
    f3 = x[5].copy()
    f4 = (x[4]*B[2] - x[5]*B[1]) * q/m
    f5 = (x[5] * B[0] - x[3] * B[2]) * q/m
    f6 = (x[3] * B[1] - x[4] * B[0]) * q/m

    return [f1,f2,f3,f4,f5,f6]

# the complete function for solving the bewegungsgleichung of a charges particle in a magnetic field
def solveSystem(init_x, init_v, b_field,dt,nr,m,q, B1_fac = 0):
    # for e)

    gradient_field = b_field.copy()

    # --start--
    times = np.arange(start=dt,stop=nr*dt+dt,step=dt)
    h = dt                                     # stepsize

    # starting vector
    X0 = [init_x[0],init_x[1],init_x[2],init_v[0],init_v[1],init_v[2]]
    X = [[0 for i in range(len(X0))] for j in range(len(times))]
    X[0] = X0

    # solve with eRK4 method
    A3 = [[0,0,0,0],[0.5,0,0,0],[0,0.5,0,0],[0,0,1,0]]
    b3 = [1/6,1/3,1/3,1/6]
    c3 = [0,1/2,1/2,1]

    # calculate Xn(t)
    for n in range(len(X) - 1):
        # arange the ki vector
        ki = [[0 for i in range(len(X0))] for i in range(4)]

        # the gradient field ( gets stronger with distance)
        d = float(np.linalg.norm(np.array([X[n][0],X[n][1], X[n][2]]) - np.linalg.norm(init_x)))
        gradient_field[0] = b_field[0] * (1+B1_fac*d)
        gradient_field[1] = b_field[1] * (1 + B1_fac * d)
        gradient_field[2] = b_field[2] * (1 + B1_fac * d)

        for i in range(4):
            t_arg = times[n] + c3[i] * h
            x_arg = X[n] + h * np.dot(A3[i], ki)
            ki[i] = func1(x_arg, t_arg, m, q, gradient_field)

        X[n + 1] = (X[n] + h * np.dot(b3, ki))

    return X

# estimates gyro radius and gyro frequency, i: index [0..3] , values from simulation
def estimateRadiusAndFrequency(X,i):
    # estimate mid-point of helix (average)
    xyz_vec = np.transpose(X)[0:3]          # x1 ... x3 are positions
    vx_vy_vz_vec = np.transpose(X)[3:6]       # x4..x6 are velocities

    [x_mid, y_mid, z_mid] = [np.average(xyz_vec[0]), np.average(xyz_vec[1]),
                             np.average(xyz_vec[2])]  # avg position of simulation

    [vx, vy, vz] = [np.average(vx_vy_vz_vec[0]), np.average(vx_vy_vz_vec[1]),
                    np.average(vx_vy_vz_vec[2])]  # avg speeds of simulation

    v_parallel = np.linalg.norm([vx, vy, vz])  # average parallel speed
    # get median of calculated helix
    xyz_vec = np.transpose(X)[0:3]
    [x_med, y_med, z_med] = [xyz_vec[0][int(_nr_steps / 2)], xyz_vec[1][int(_nr_steps / 2)],
                             xyz_vec[2][int(_nr_steps / 2)]]

    # subtract both vectors for radius
    radius_vektor = [x_mid - x_med, y_mid - y_med, z_mid - z_med]
    r_estimate = np.linalg.norm(radius_vektor)  # absolute size of radius vektor

    v_tang = np.sqrt(_V0s[i] ** 2 - v_parallel ** 2)  # tangential part (pythagoras)
    w_estimate = v_tang / r_estimate

    return r_estimate,w_estimate

# method to calculate gyro radius and gyro frequency, i: index [0..3]
def calcRadiusAndFrequency(i):
    v_parallel = np.dot(_ev, _eb) * _V0s[i]  # calc parallel-part of velocity ( along B)
    v_tang = np.sqrt(_V0s[i] ** 2 - v_parallel ** 2)  # tangential part (pythagoras)
    r_g = v_tang * _m1 / (_q1 * _B0s[i])  # radius = V_t*m /qB
    w_g = v_tang/r_g
    return r_g,w_g

# ----------- START ----------------------------------------------------------

# --- a) --- solve system using the 3 methods euler, midpoint, eRK4
print("\n --------------------------- a) ------------------------ ")

x0 = [0,0,0]        # starting position
_m1 = 1.6735e-27     # mass of particle
_q1 = 1.6022e-19     # charge of particle

_dt_factor = 0.001   # factor by which Tau is multiplied for dt
_nr_steps = 2e4      # number of steps to calculate

_eb = [np.sqrt(3) / 4, 1 / 4, np.sqrt(3) / 2]      #direction of B
_ev = [np.sqrt(3) / 2, 1 / 2, 0]     # direction of V


_B0s = [3, 0.1, 2.5e-8, 1e-9]    # B fields in [T]
E_KIN = [1e5,0.3,50,0.016]    # kinetic energy in [eV]

# --- b) --- estimate optimal timestep
print("\n --------------------------- b) ------------------------ ")

_V0s = E_KIN.copy()
Tau = E_KIN.copy()          # calculate 1 Tau of revolve
for i in range(len(_V0s)):
    _V0s[i] = np.sqrt(2 * E_KIN[i] * _q1 / _m1)
    Tau[i] = 2 * np.pi * _m1 / (_q1 * _B0s[i])

dt = np.multiply(Tau, _dt_factor)

print("Tau [s]: ",Tau)
print("v0 [m/s]: ", _V0s)
print("dt [s]: ", dt)


# --- c) --- start calculation with B , v , ..
print("\n --------------------------- c) ------------------------ ")

B0_matrix = [[ 0 for i in range(3)] for j in range(len(_B0s))]
V0_matrix = B0_matrix.copy()
for i in range(len(B0_matrix)):
    B0_matrix[i] = np.multiply(_eb, _B0s[i])
    V0_matrix[i] = np.multiply(_ev, _V0s[i])

# obtain solutions and plot
Xs = []
titles = ["Fusion Device", "Sunspot", "Earth magnetosphere", "interplanetary medium"]
for i in range(len(B0_matrix)):
    Xs.append(solveSystem(x0, V0_matrix[i], B0_matrix[i], dt[i], _nr_steps, _m1, _q1))

for i in range(len(B0_matrix)):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(np.transpose(Xs[i])[0], np.transpose(Xs[i])[1], np.transpose(Xs[i])[2])
    plt.xlabel("X / m")
    plt.ylabel("Y / m")
    plt.title(titles[i])
    plt.show()

# ---- d) ----- from simulation: estimate values of gyroradius and gyrofrequency and compare
print("\n --------------------------- d) ------------------------ ")

for i in range(len(Xs)):
    re,we = estimateRadiusAndFrequency(Xs[i],i)
    rc,wc = calcRadiusAndFrequency(i)
    print("i:",i,"\nr_g_estimate [m]: ",re,"\tw_g_estimate: [rad/s] ",we,"\nr_g_calc: ",rc,"\tw_g_calc",wc)
    print("\n")

# ---- e) -----  B field with gradient
print("\n --------------------------- e) ------------------------ ")
B1_factor = 1   # effekt sehr gut sichtbar

# fusion device with gradient
X_gradient = solveSystem(x0, V0_matrix[0], B0_matrix[0], dt[0], _nr_steps, _m1, _q1,B1_factor)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(np.transpose(X_gradient)[0], np.transpose(X_gradient)[1], np.transpose(X_gradient)[2])
ax.plot(np.transpose(Xs[0])[0], np.transpose(Xs[0])[1], np.transpose(Xs[0])[2])
plt.xlabel("X / m")
plt.ylabel("Y / m")
plt.legend(('X with gradient-field','Xs'))
plt.title("particle in Field with gradient")
plt.show()
