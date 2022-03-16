import random

from datetime import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# TODO calculate electronic band strucutre of some material using the plane wave method / tight binding method

n = 1e28   # electronic density





kb = 1.38e-23
h = 6.626e-34
u = 1.66e-27


def S(N,V,U,m):
    return N*kb*(np.log(V/N * (4*np.pi*m*U/(3*N*h**2))**(1.5))+5/2)

print("SA = ",S(2.17e+26,3,4.5e+5,40*u))
print("SC = ",S(2.96e+26,4,1.5*kb*1000*2.96e+26,20*u))