import numpy as np
import matplotlib.pyplot as plt

# DESCRIPTION
# Intrinsic Semiconductor: calculate hole concentration

#Global Variables

_kb = 1.380649e-23      #J/K
_T0 = 300               #K
_hbar = 1.05457e-34     #J*s
_e = 1.602e-19          #A*s
_ev_to_J = 1.60218e-19  #J / ev
_me = 9.109e-31         #kg .. mass of electron
_PI = np.pi

# Conductivity Tensor
def calc_sigma(n,mu_n, p, mu_p):
    return n*_e*mu_n + p*_e*mu_p;

# Mobility Constant
def calc_mu_p(v_drift, E_ext):
    return v_drift/E_ext;

# Diffusion const
def calc_Dp(mup,T):
    return _kb*T * mup /_e;


# --------Calculate Sigma------------
ni = 1.5e+10 # 1/cm3
mup = 450 #cm2/Vs
mun = 1500 # cm2 / Vs

Na = 2e+17 #1/cm3
p1 = Na

n1 = ni**2 / p1  #1/cm3

print(f'n1: {n1:e}')
print("Sigma in 1/Ohm*cm is ",calc_sigma(n1,mun,p1,mup))

# ----------- Calculate Drift mobility constant

s0 = 0.01 #m
t0 = 0.0001 #s
v_drift = s0/t0 #m/s
E0 = 7311 # V / m

mup0 = calc_mu_p(v_drift,E0) # m2/Vs
print("mobility const cm2/Vs is ",mup0*1e+4)
print("diffusion const cm2/s is ",calc_Dp(mup0*1e+4,300))
