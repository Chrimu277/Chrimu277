# Muehleder Christoph 01604413

import numpy as np
import matplotlib.pyplot as plt

# --------calculate functions for exercise-------

_X1 = np.e      #nulldurchgang bei X1, X2
_X2 = np.pi

# für y1(x), berechne phi, omega
_OMEGA = _X2 / (_X2 - _X1)
_PHI = -_X1 * _X2 / (_X2 - _X1)

# für y2(x) , berechne x0, d
_D = 10      # parabelkuppel
_X0 = (_X1 + _X2) / 2
_A = _D / (_X1 / 2 - _X2 / 2) ** 2

xs = np.arange(0,5,0.01)
y1 = np.sin(xs * _OMEGA + _PHI)
y2 = _A * (xs - _X0) ** 2 - _D

plt.figure()
plt.plot(xs, y2, 'r', label = 'square') # plotting t, a separately
plt.plot(xs, y1, 'b', label = 'sine') # plotting t, b separately
plt.ylim(-_D, _D)
plt.xlim(0,np.max(xs))
plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("y1(x), y2(x)")
plt.show()




# ------FUNCTIONS DEFINITIONS--------
# estimate value of pi and e for y1(x) , y2(x)

# y1(x) , y1(pi) = 0 , y1(e) = 0 , sine function
def y1(x,w,phi):
    return np.sin(w*x+phi);

# first derivate of y1 by x
def y1_der1(x,w,phi):
    return w*np.cos(w*x+phi)

# y2(x) , y2(pi) = 0 , y2(e) = 0 , square function
def y2(x, a,x0,d):
    return a*(x-x0)**2 - d

def y2_der1(x, a,x0,d):
    return 2*a*(x-x0);

# newton-raphson method for finding zeroes of y1(x) , 1 starting value
def getZero_newton_raphson1(x_start, maxsteps = 10, digits = 6):
    xs = []
    for i in range(maxsteps):
        x_start = x_start - y1(x_start, _OMEGA, _PHI) / y1_der1(x_start, _OMEGA, _PHI)
        xs.append(x_start)
        if i > 1 and check_equality(xs[i],xs[i-1], digits):     # check if numbers have same [digits] after comma digits
            break;
    return xs

# newton-raphson method for finding zeroes of y2(x) , 1 starting value
def getZero_newton_raphson2(x_start, maxsteps = 10, digits = 6):
    xs = []
    for i in range(maxsteps):
        x_start = x_start - y2(x_start,_A,_X0,_D) / y2_der1(x_start,_A,_X0,_D)
        xs.append(x_start)
        if i > 1 and check_equality(xs[i],xs[i-1], digits):     # check if numbers have same [digits] after comma digits
            break;
    return xs

# regula falsi method, xs contains 2 start values, for y1(x)
def getZero_RegulaFalsi1(xs, maxsteps = 10, digits = 6):
    for i in range(1, maxsteps - 1):
        if check_equality(xs[i],xs[i-1], digits):               # check if numbers have same [digits] after comma digits
            break;
        xnew = xs[i] - (y1(xs[i],_OMEGA, _PHI) * (xs[i] - xs[i - 1]) / ((y1(xs[i], _OMEGA, _PHI) - y1(xs[i - 1], _OMEGA, _PHI))))
        xs.append(xnew)
    return xs

# regula falsi method, xs contains 2 start values, for y2(x)
def getZero_RegulaFalsi2(xs, maxsteps = 10, digits = 6):
    for i in range(1,maxsteps-1):
        if check_equality(xs[i],xs[i-1], digits):
            break;
        xnew = xs[i] - (y2(xs[i],_A,_X0,_D) * (xs[i] - xs[i-1]) / ((y2(xs[i],_A,_X0,_D) - y2(xs[i-1],_A,_X0,_D))))
        xs.append(xnew)
    return xs

# bisection method with interval [a0,b0] for y1(x)
def getZero_Bisection1(a0,b0, maxsteps = 10, digits = 6):
    xs = []
    for i in range(maxsteps):
        x0 = (a0+b0) / 2
        xs.append(x0)

        if (y1(x0,_OMEGA,_PHI) * y1(a0,_OMEGA,_PHI) < 0):
            b0 = x0
        elif (y1(x0,_OMEGA,_PHI) * y1(a0,_OMEGA,_PHI) > 0) :
            a0 = x0;
        if i > 1 and check_equality(xs[i],xs[i-1], digits):
            break;
    return xs

# bisection method with interval [a0,b0] for y2(x)
def getZero_Bisection2(a0, b0, maxsteps=10, digits = 6):
    xs = []
    for i in range(maxsteps):
        x0 = (a0 + b0) / 2
        xs.append(x0)

        if (y2(x0, _A, _X0, _D) * y2(a0, _A, _X0, _D) < 0):
            b0 = x0
        elif (y2(x0, _A, _X0, _D) * y2(a0, _A, _X0, _D) > 0):
            a0 = x0;
        if i > 1 and check_equality(xs[i],xs[i-1], digits):
            break;
    return xs

# get division of 1/A without divisions!
def getDivision_1_A(x_start,a, maxsteps = 10, digits = 6):
    xs = []
    xs.append(x_start)
    for i in range(maxsteps):
        # newton raphson method on f(x)
        x_start = 2*x_start - a*(x_start**2)
        xs.append(x_start)

        if i > 1 and check_equality(xs[i], xs[i - 1], digits):
            break;
    return xs

# this method checks if 2 floats are the same until commadigit [digits]
def check_equality(x1,x2,comma_digits):
    if (int)(x1*np.power(10,comma_digits) - x2*np.power(10,comma_digits)) == 0:
        return True;



# -------START PROGRAM --------------

_COMMA_DIGITS = 6
_MAXSTEPS = 10

#estimate E with newton-raphson method for both functions
e_1_nr = getZero_newton_raphson1(2.6,_MAXSTEPS,_COMMA_DIGITS) # max 100 iterations
e_2_nr = getZero_newton_raphson2(2.8,_MAXSTEPS,_COMMA_DIGITS)

#estimate PI with newton-raphson method for both functions
pi_1_nr = getZero_newton_raphson1(3.1,_MAXSTEPS,_COMMA_DIGITS)
pi_2_nr = getZero_newton_raphson2(3.2,_MAXSTEPS,_COMMA_DIGITS)

#estimate E with regula falsi method for both functions
e_1_rf = getZero_RegulaFalsi1([2.6,2.8],_MAXSTEPS,_COMMA_DIGITS);
e_2_rf = getZero_RegulaFalsi2([2.6,2.8],_MAXSTEPS,_COMMA_DIGITS);

#estimate PI with regula falsi for both functions
pi_1_rf = getZero_RegulaFalsi1([3.1,3.2],_MAXSTEPS,_COMMA_DIGITS);
pi_2_rf = getZero_RegulaFalsi2([3.1,3.2],_MAXSTEPS,_COMMA_DIGITS);

#estimate E with bisection method for both functions
e_1_bs = getZero_Bisection1(2.6,2.8,_MAXSTEPS,_COMMA_DIGITS)
e_2_bs = getZero_Bisection2(2.6,2.8,_MAXSTEPS,_COMMA_DIGITS)

#estimate PI with bisection method for both functions
pi_1_bs = getZero_Bisection1(3.1,3.2,_MAXSTEPS,_COMMA_DIGITS)
pi_2_bs = getZero_Bisection2(3.1,3.2,_MAXSTEPS,_COMMA_DIGITS)


print(e_1_nr)
print(e_2_nr)
print(pi_1_nr)
print(pi_2_nr)

print(e_1_rf)
print(e_2_rf)
print(pi_1_rf)
print(pi_2_rf)

print(e_1_bs)
print(e_2_bs)

print(pi_1_bs)
print(pi_2_bs)


plt.figure()
plt.plot(range(len(e_1_nr)),e_1_nr, label = 'e from y1, NR, x0 = 2.6')
plt.plot(range(len(e_2_nr)),e_2_nr, label = 'e from y2, NR, x0 = 2.8')
plt.plot(range(len(e_1_rf)),e_1_rf, label = 'e from y1, RF, x0,1 = 2.6,2.8')
plt.plot(range(len(e_2_rf)),e_2_rf, label = 'e from y2, RF, x0,1 = 2.6,2.8')
plt.plot(range(len(e_1_bs)),e_1_bs, label = 'e from y1, BS, x0,1 = 2.6,2.8')
plt.plot(range(len(e_2_bs)),e_2_bs, label = 'e from y2, BS, x0,1 = 2.6,2.8')
plt.legend()
plt.show()

plt.plot(range(len(pi_1_nr)),pi_1_nr, label = 'e from y1, NR, x0 = 3.1')
plt.plot(range(len(pi_2_nr)),pi_2_nr, label = 'e from y2, NR, x0 = 3.2')
plt.plot(range(len(pi_1_rf)),pi_1_rf, label = 'e from y1, RF, x0,1 = 3.1,3.2')
plt.plot(range(len(pi_2_rf)),pi_2_rf, label = 'e from y2, RF, x0,1 = 3.1,3.2')
plt.plot(range(len(pi_1_bs)),pi_1_bs, label = 'e from y1, BS, x0,1 = 3.1,3.2')
plt.plot(range(len(pi_2_bs)),pi_2_bs, label = 'e from y2, BS, x0,1 = 3.1,3.2')
plt.legend()
plt.show()

# g) calc division without using divisions, make algorithm with NR method
print("division")
print(getDivision_1_A(0.1,7))       # calculate 1/7 starting at 0.1

#plot iteration process