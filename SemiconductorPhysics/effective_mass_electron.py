import numpy as np
import matplotlib.pyplot as plt

# DESCRIPTION

#Global Variables

_kb = 1.380649e-23      #J/K
_T0 = 300               #K
_hbar = 1.05457e-34     #J*s
_e = 1.602e-19          #A*s
_ev_to_J = 1.60218e-19  #J / ev
_me = 9.109e-31         #kg .. mass of electron


m_dash = 0.2*_me # eletron masses
k1 = 6e+8           # 1/m

E = (_hbar*k1)**2/(2*m_dash)
E_ev = E/_ev_to_J

print(E_ev)