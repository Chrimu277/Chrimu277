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
_m_He = 6.6e-27

#Electron Energy eV -> Lambda nm
# E = p**2 / 2m = hbar2 k2 / 2m
# p = hbar k = h / lambda
def get_lambda_kineticmass(energy,mass):
    energy *= _ev_to_J
    lam = 2*np.pi*_hbar / np.sqrt(2*mass*energy)  # 1 Option
    #lam = np.sqrt(_hbar**2 * np.pi**2 * 2 /(_me*energy)) #2 Option
    return lam;

# Lambda --> Energy
# formulas same
def get_energy_kineticmass(lam,mass):
    energy = (2*np.pi/lam)**2 * _hbar**2 / (2*mass)
    return energy;

#Energy ev -> Lambda m
def get_lambda_phonton(energy):
    lam = 1240 / energy
    #lam = np.sqrt(_hbar**2 * np.pi**2 * 2 /(_me*energy)) #2 Option
    return lam*1e-9;

# Lambda nm --> Energy ev
# formulas same
def get_energy_phonton(lam):
    energy = (2*np.pi/lam)**2 * _hbar**2 / (2*_me)
    return energy;

# BSP:
E1 = 0.2 #eV

lambda_electron  = get_lambda_kineticmass(E1,_me)
lambda_helium = get_lambda_kineticmass(E1, _m_He)
lambda_phonon = get_lambda_phonton(E1)

print("l el:",lambda_electron,"\nl He: ",lambda_helium, "\nl Pho:0",lambda_phonon)