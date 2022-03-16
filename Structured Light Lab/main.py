import numpy as np
import matplotlib.pyplot as plt


def In(t,A,B,C,D,t0):
    theta = (t+t0)*np.pi/180
    return 0.5*(A+B*np.sin(2*theta) + C*np.cos(2*theta) + D*np.sin(4*theta))

def A(I, theta):
    if len(I) != len(theta):
        print("err")
        return -1
    A = 0.0
    for i in range(len(I)):
        A += I[i]
    return A * 2/len(I)

def B(I, theta):

    B = 0.0

    if len(I) != len(theta):
        print("err")
        return -1

    for i in range(len(I)):
        B += I[i]*np.sin(2*theta[i]*np.pi/180)

    return B * 4/len(I)

def C(I, theta):
    C = 0.0
    if len(I) != len(theta):
        print("err")
        return -1

    for i in range(len(I)):
        C += I[i] * np.cos(4 * theta[i]*np.pi/180)
    return C * 4/len(I)

def D(I, theta):

    D = 0.0

    if len(I) != len(theta):
        print("err")
        return -1

    for i in range(len(I)):
        D += I[i] * np.sin(4 * theta[i]*np.pi/180)

    return D * 4/len(I)

def stokesvektor(a,b,c,d):
    sv = [a-c,2*c,2*d,b]
    sv = sv / sv[0]
    return sv

def norm(vec):
    norm = 0.0
    for i in range(len(vec)):
        norm += vec[i]**2

    return np.sqrt(norm)

# theta data
theta1 = np.arange(0.0,180.0,10)
theta2 = np.arange(0.0,180.0,10)
theta0_1 = 9
theta0_2 = 79

# intensity data
I1 = np.array([17.4, # linear analysis
      39.6,
      66.9,
      86.5,
      87.0,
      70.8,
      42.2,
      18.9,
      8.1,
      16.1,
      41.5,
      68.9,
      86.6,
      87.8,
      68.5,
      44.3,
      18.2,
      8.1])
I2 = np.array([85.5, # circular analysis
      114.3,
      138.6,
      151.2,
      151.2,
      138.8,
      114.3,
      89.0,
      63.7,
      42.1,
      27.3,
      18.4,
      13.6,
      12.1,
      14.3,
      21.8,
      35.7,
      58.5])
error = 0.1 # 10% up
I1 = I1 * (1+error)
I2 = I2 * (1+error)

Theta0_I10 = theta1[np.nanargmin(I1)]
Theta0_I20 = theta2[np.nanargmin(I2)]
print("theta 0#s could be: ",Theta0_I10, Theta0_I20)

N1 = len(I1)
N2 = len(I2)

theta1 += theta0_1
theta2 += theta0_2

A1 = A(I1,theta1)
A2 = A(I2, theta2)

B1 = B(I1,theta1)
B2 = B(I2, theta2)

C1 = C(I1,theta1)
C2 = C(I2, theta2)

D1 = D(I1,theta1)
D2 = D(I2, theta2)

print("CALCULATION BY FORMULA")
print("ABCD_lin = ",A1,B1,C1,D1," with theta0 = ",theta0_1)
print("ABCD_lin = ",A2,B2,C2,D2," with theta0 = ",theta0_2)
s1 = stokesvektor(A1,B1,C1,D1)
s2 = stokesvektor(A2,B2,C2,D2)
print("S_lin = ",s1)
print("S_circ = ",s2)
print("P1 = ",norm(s1[1:])/s1[0])
print("P2 = ",norm(s2[1:])/s2[0])


# CALCULATE PER FIT
print("CALCULATION BY FIT")
from scipy.optimize import curve_fit
popt, pcov = curve_fit(In, np.arange(0,180,10), I1,p0=[A1,B1,C1,D1,Theta0_I10],check_finite=True)
popt2, pcov2 = curve_fit(In, np.arange(0,180,10), I2,p0=[A2,B2,C2,D2,Theta0_I20],check_finite=True)

fit_In_lin = In(np.arange(0,180,2), *popt)
fit_In_circ = In(np.arange(0,180,2), *popt2)

err1 = 0.1 * I1
err2 = 0.1 * I2

plt.plot(np.arange(0,180,2)+theta0_1,fit_In_lin)
plt.errorbar(theta1, I1,c='black',yerr= err1, ecolor='red',fmt='o')
plt.title("Fit of Intensity Measurement of Linear Polarization")
plt.xlabel("Theta / °")
plt.ylabel("I / uW")
plt.legend(["Fit curve","measured values"])
plt.grid()
plt.show()

plt.plot(np.arange(0,180,2)+theta0_2,fit_In_circ)
plt.errorbar(theta2, I2,c='black',yerr= err2, ecolor='red',fmt='o')
plt.title("Fit of Intensity Measurement of Circular Polarization")
plt.xlabel("Theta / °")
plt.ylabel("I / uW")
plt.legend(["Fit curve","measured values"])
plt.grid()
plt.show()

print(len(popt),"popt_lin (A,B,C,D, t0): ",*popt)
print(len(popt2),"popt_lin (A,B,C,D, t0): ",*popt2)

print("S_lin = ",stokesvektor(*popt[:4]))
print("S_circ = ",stokesvektor(*popt2[:4]))

