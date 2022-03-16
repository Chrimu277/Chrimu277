# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------


# runge kutta numerical method to solve differential equation
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
    #TODO untere dreiecksmatrix

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



# given function on problem a)
_alpha = -0.1
def func1(y,t):
    r = np.multiply(_alpha,y)
    return r


# analytical solution of yd(t)
def yd_analytical(y0,t,a):
    return y0*np.exp(a*t)


# solved function on problem f)
def func2(y,t):
    f = np.array([0,0])
    f[0] = y[1]
    f[1] = 4*t - 2 * y[0] - 3 * y[1]
    return f;


# analytical solution of yf(t)
def yf_analytical(t):
    return 2*t - 3 - 4*np.exp(-2*t) + 8*np.exp(-1*t)

# ----------- START ----------------------------------------------------------

# --- a) --- runge kutta method: solve linear equation system---------------------------------

# --- b) ---- check if butcher tableau is correct ----

# --- c,d) ---------------------------------------------------
print("\n --------------------------- c) d) -----------------------------")
#y' = a * y(t)

# A ,b ,c for Euler-method solving
A = [[0]]
b = np.array([1])
c = np.array([0])

# initial values
_Y0 = np.array([2])

# time array
dt = 1
dt2 = 0.01
tmax = 50
times = np.arange(dt2,tmax+dt,dt)
times2 = np.arange(dt2,tmax+dt2,dt2)

# calculate with dt =1s , dt2 = 0.01s
yd1 = RungeKutta(func1,_Y0,A,b,c,times)
yd2 = RungeKutta(func1,_Y0,A,b,c,times2)


# ---- e)  solve analytically and plot alltogether ------------------
print("\n --------------------------- e) -----------------------------")
yd3 = yd_analytical(_Y0,times2,_alpha)

#plot
plt.figure()
plt.plot(times,yd1)
plt.plot(times2,yd2)
plt.plot(times2,yd3)
plt.legend(('dt = 1',' dt = 0.01' , 'analytical solution'))
plt.title("Runge Kutta: Euler Method (order 1)")
plt.grid()
plt.show()


# ---- f,g ) --------- solve differential equation of 2nd order f(2)  dt = 0.01 , tmax = 20
print("\n --------------------------- f) g) -----------------------------")

# inital values of y(t), y'(t)
_Y0 = [1,2]

# time
dt = 0.001
tmax = 20
times = np.arange(0,tmax,dt)

# solve with euler
yf1 = RungeKutta(func2,_Y0,A,b,c,times)

#solve with explicit midpoint method
A2 = [[0,0],[0.5,0]]
b2 = [1/2,1/2]
c2 = [0,1/2]
yf2 = RungeKutta(func2,_Y0,A2,b2,c2,times)

# solve with eRK4 method
A3 = [[0,0,0,0],[0.5,0,0,0],[0,0.5,0,0],[0,0,1,0]]
b3 = [1/6,1/3,1/3,1/6]
c3 = [0,1/2,1/2,1]
yf3 = RungeKutta(func2,_Y0,A3,b3,c3,times)

# get solution from hand-calculated solution
yf4 = yf_analytical(times)

# plot alltogether
plt.figure()
plt.plot(times,yf1)
plt.plot(times,yf2)
plt.plot(times,yf3)
plt.plot(times,yf4)
plt.legend(('yf1(t) Euler','d/dt yf1(t)' ,'yf2(t) eMidpoint','d/dt yf2(t)' ,'yf3(t) eRK4','d/dt yf3(t)', 'yf4(t) Analyt.'))
plt.title("Runge Kutta: Euler Method, E-Midpoint, E-RK4, Analytical Solution")
plt.grid()
plt.show()


