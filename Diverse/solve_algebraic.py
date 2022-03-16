# SOURCES :
#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

#solve equation for T
# R = (1-e^-T) / (1-e^(-aT))
# or rightside - leftside = 0!
a = 0.5
R = 1.6

#define lambda function
func = lambda T : R - ((1.0 - np.exp(-T))/(1.0 - np.exp(-a*T)))
T_arr = np.linspace(-5,5,1000)

# plot expression
plt.plot(T_arr,func(T_arr))
plt.xlabel("Tau")
plt.ylabel("expression value -> must go to 0")
plt.grid()
plt.show()

#solve and test
T0_guess = 2
T_solve = fsolve(func, T0_guess)
print("type of solved variable:",type(T_solve))
print("Solution found for Tau",T_solve,sep=" = ")
print("f(",T_solve[0],") =",func(T_solve[0]))

# ------ END ----------------