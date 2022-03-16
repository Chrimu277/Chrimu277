# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------

# runge kutta method for numerically solving differential equations
def RungeKutta(function, y0, A, b, c, time_vec):

    y = [[0 for i in range(len(y0))] for j in range(len(time_vec))]
    y[0] = y0

    h = time_vec[2]- time_vec[1]        #stepsize

    valid,order = checkButcherTableau(A,b,c)    #check A,b,c
    if not valid:                               # bye bye
        print("butcher not valid")
        return np.zeros(len(time_vec));

    # valid
    print("order is ",order)

    if order > 1:
        #y(n+1) = y(n) + h * sum_i bi ki
        for n in range(len(y)-1):
            # arange the ki vector
            ki = [[0 for  i in range(len(y0))] for i in range(order)]
            for i in range(order):
                t_arg = time_vec[n]+c[i]*h
                y_arg = y[n] + h* np.dot(A[i],ki)
                ki[i] = (function(y_arg,t_arg))

            y[n+1] = (y[n] + h* np.dot(b,ki))

    # order = 1, euler method
    else:
        for n in range(len(y)-1):
            y[n+1] = y[n] + h * function(y[n],time_vec[n])

    return y


# checks if butcher tableau is correct and valid
def checkButcherTableau(A,b,c):

    valid = True;
    tol = 1e-9;     #tolerance because 1 ~ 0.99999999999999

    #first element of c must be 0
    if c[0] != 0:
        print("Error: c0 not 0!")
        valid = False

    # sum of b must be 1
    if np.abs(np.sum(b)-1) > tol:
        print("Error: sum of b  not 1!")
        print(np.sum(b))
        valid = False

    # the dimensions must fit
    if len(A) != len(A) or len(c) != len(b) or len(c) != len(A):
        print("Dimensions do not fit")
        print("A: ",np.shape(A),", b: ",np.shape(b), ", c: ",np.shape(c))
        valid = False;

    # each sum of row i of A must be ci
    for i in range(len(A)):
        if np.abs(np.sum(A[i])-c[i]) > tol:
            print("sum A",i," not c",i)
            print(np.sum(A[i])," and ", c[i])
            valid = False;

    return valid, len(b)


# given function for with Y' = f(Y,t..)
def func1(y,t):
    f1 = y[1]
    f2 = -(_G/_L)*np.sin(y[0])
    return np.array([f1,f2])


# calculate the hamiltonian (energy) of the problem
def Hamiltonian(y): #y[0] = phi , y[1] = w
    return _M*((y[1]*_L)**2)/2 - _M*_G*_L * np.cos(y[0])


# plots the data with title and legends
def plotSolution(data,title,legend):
    plt.figure()
    for i in range(len(data)):
        plt.plot(np.transpose(data[i])[0], np.transpose(data[i])[1])
    plt.title(title)
    plt.xlabel("phi / rad")
    plt.ylabel("omega / rad/s")
    plt.legend(legend)
    plt.grid()
    plt.show()
    return
# ----------- START ----------------------------------------------------------

# --- a) --- solve system using the 3 methods euler, midpoint, eRK4
print("\n ----- a) ------")

_G = 9.81               # earths acceleartion
_L = 1                  # length of pendulum
_PHI0 = np.pi / 3       # initial phi
_M = 1                  # mass of pendulum mass

_Y0 = [_PHI0,0]         # initial vektor

# time
dt = 0.01
tmax = 5
times = np.arange(dt,tmax+dt,dt)

# A ,b ,c for Euler-method solving
A1 = [[0]]
b1 = np.array([1])
c1 = np.array([0])
ya1 = RungeKutta(func1,_Y0,A1,b1,c1,times)  # contains phi, omega

# solve with explicit midpoint method
A2 = [[0,0],[0.5,0]]
b2 = [1/2,1/2]
c2 = [0,1/2]
ya2 = RungeKutta(func1,_Y0,A2,b2,c2,times)

# solve with eRK4 method
A3 = [[0,0,0,0],[0.5,0,0,0],[0,0.5,0,0],[0,0,1,0]]
b3 = [1/6,1/3,1/3,1/6]
c3 = [0,1/2,1/2,1]
ya3 = RungeKutta(func1,_Y0,A3,b3,c3,times)

# --- b) ---- plot trajectories over time
print("\n ----- b) ------")
plt.figure()
plt.plot(times,ya1)
plt.plot(times,ya2)
plt.plot(times,ya3)
plt.legend(('phi1(t) Euler','w1(t)' ,'phi2(t) eMidpoint','w2(t)' ,'phi3(t) eRK4','w(t)'))
plt.title("Pendulum: Euler Method, E-Midpoint, E-RK4")
plt.xlabel("time / s")
plt.ylabel("phi(t), w(t) / s")
plt.grid()
plt.show()
plt.show()


# --- c) ---- plot hamiltonian (total energy of system)
print("\n ----- c) ------")

h1 = Hamiltonian(np.transpose(ya1))
h2 = Hamiltonian(np.transpose(ya2))
h3 = Hamiltonian(np.transpose(ya3))

plt.figure()
plt.plot(times,h1)
plt.plot(times,h2)
plt.plot(times,h3)
plt.legend(('H1(t)','H2(t)','H3(t)'))
plt.title("Plot of Hamiltonians")
plt.xlabel("time / s")
plt.ylabel("E_TOT(t) [J]")
plt.grid()
plt.show()

# ---- d) ---- plot system dynamics phi(t) against w(t)
print("\n ----- d) ------")
# 10 equally distributed omegas, solve system for each
_W0 = np.linspace(start= -np.sqrt(2*_G/_L) , stop = np.sqrt(2*_G/_L) , num=10)
_Y0 = np.transpose([[_PHI0 for i in range(len(_W0))],_W0])         # initial vektors

#get all the solutions
Ys_erk4 = []
Ys_emp = []
Ys_euler = []
legends = []

for i in range(len(_W0)):
    Ys_euler.append(RungeKutta(func1,_Y0[i],A1,b1,c1,times))
    Ys_emp.append(RungeKutta(func1, _Y0[i], A2, b2, c2, times))
    Ys_erk4.append(RungeKutta(func1, _Y0[i], A3, b3, c3, times))
    legends.append("w0 = "+str(_Y0[i][1])+", phi0 = "+str(_Y0[i][0]))

# plot all solutions
plotSolution(Ys_erk4,"Plot of System Dynamics, eRK4",legends)
plotSolution(Ys_emp,"Plot of System Dynamics, Euler Method",legends)
plotSolution(Ys_erk4,"Plot of System Dynamics, eRK4",legends)

