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

#function definitions

#fermi dirac probability
def fermi_dirac(E,T,Ef = 0):
    t1 = np.exp((E-Ef)/(_kb*T))
    p = 1 / (1 + t1)
    return p


#runtime variables

#make arrays for temperatures, Energies, probabilities
T_arr = np.arange(start=100, stop=2500, step = 100)
E_arr = np.linspace(-1,1,1000)
E_arr *= _ev_to_J
P_arr = []

# calculate probability of occupation
for T in T_arr:
    P_arr.append(fermi_dirac(E_arr ,T))

# plot all
for P in P_arr:
    plt.plot(E_arr, P)

plt.xlabel("E / J")
plt.ylabel("P_occ(E,T)")
plt.grid(which='both')
plt.show()


#calculate 1 specific probability
T = 206 #K
E = -0.07*_ev_to_J #J
Ef = 0;

P_occ = fermi_dirac(E,T)
print(f'Probability of occupation at T = {T} ; Ef = {Ef}; E = {E} \n is: {P_occ}')