import random
from datetime import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# ideales Gas:
# pv = nRT // n .. number of moles   R = kB NA .. gas constant    T .. Temperature
# Ekin = f/2 kb T  f.. DOF

#nature constants ---
kb = 1.38e-23 #J/K
NA = 6.022e+23 # 1, avogadros number
R = kb*NA # J/K
u = 1.66e-27 #kg, atomic mass unit
# ----

# -----------------------------------------
# calculate gas internal energy (kinetic)
# return unit: Joule
# n: nr of moles
# T: temperature
# f: DOF
def E_tot(f,n,T):
    return f/2*kb*T*n*NA
# -----------------------------------------


if __name__ == '__main__':
    print("START")
    print("E_A:",E_tot(f=3,n= 3.608e-3, T= 100))
    print("E_B:",E_tot(f=6,n= 1e-5, T = 1200))