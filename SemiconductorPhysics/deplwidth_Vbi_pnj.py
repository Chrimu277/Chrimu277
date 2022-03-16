import numpy as np
import matplotlib.pyplot as plt

# Christoph Mühleder

# nature constants
_kb = 1.380649e-23      #J/K
_T0 = 300               #K
_hbar = 1.05457e-34     #J*s
_e = 1.602e-19          #A*s
_ev_to_J = 1.60218e-19  #J / ev
_me = 9.109e-31         #kg .. mass of electron
_eps0 = 8.8541878128e-12 #F/m  vacuum permittivity

#function definitions

# calculate Depletion Width W of pn-junction
def W(eps, Nd, Na, Vbi, V):
    up = 2*eps*(Nd+Na)*(Vbi-V)
    down = _e*Nd*Na
    return np.sqrt(up/down)

def Wn_Wp(W, Na ,Nd):
    return W*Na/(Na+Nd), W*Nd/(Na+Nd)

# calculate V-built-in of pn-junction
def builtin_voltage(Nd, Na, ni, T):
    return _kb*T/_e * np.log(Nd*Na / ni**2)



# -- V builtin, Wn, Wp = ?  Voltage applied = 0 -------------
Nd = 7000000000000000 #1/cm3
Na = 9000000000000000000 #1/cm3
T1 = 300     #K
ni = 1.5e+10 #1/cm3
eps_r = 11.9

#to SI
Nd *= 1e+6
Na *= 1e+6
ni *= 1e+6

Vbi = builtin_voltage(Nd,Na,ni,T1)
print("Vbi = ",Vbi," V")
V_forward = 0

Wn, Wp = Wn_Wp(W(eps_r*_eps0,Nd,Na,Vbi,V_forward),Na,Nd)
print("Wn = ",Wn," m\nWp = ", Wp," m")