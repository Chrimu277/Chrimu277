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

# functions
def hole_conc(nc, nv, eg, t):
    p = np.sqrt(nv*nc*np.exp(-eg/(_kb*t))*((t/_T0)**(3/2)))
    return p;

#Nc, Nv in 1/m3, Eg in J, t in K
def fermi_ernergy__intr_C(nc, nv, eg, t):
    ef = eg/2 + 0.5*(_kb*t)*np.log(nv/nc)
    return ef;

# ne in 1/m3 , result Efmet in J
def fermi_energy_metal_0K(ne):
    EF_met = (_hbar ** 2 / (2 * _me)) * ((3 * np.pi ** 2 * ne) ** (2 / 3))
    return EF_met


def main():

    #runtime Variables
    Nc_300k = 2.78e+25  # 1/m3
    Nv_300k = 9.84e+24  # 1/m3
    Eg_300k = 1.12*_ev_to_J #J

    T = np.arange(100, 500)
    T2 = np.linspace(100,500)
    #print("Type of T:",type(T))
    #print("Type of T2:", type(T2))

    # pass array of T into function which wants scalar -> gives back array !!
    #p_of_T = hole_conc(Nc_300k,Nv_300k,Eg_300k,T)
    #plt.plot(T,p_of_T)
    #plt.show()

    T1 = 378
    p_T1 = hole_conc(Nc_300k,Nv_300k,Eg_300k,T1)
    print(f'hole concentration p at {T1}K is {p_T1:e}') # Todo noch falsches Ergebnis !!!!

#----------------------------------------------------

    # calculate fermi energy of Metal at T = 0K
    #WORKS
    ne = 2e+28 #electron concentration
    print("EF_met = ",fermi_energy_metal_0K(ne)," J = ", fermi_energy_metal_0K(ne)/_ev_to_J," eV")

#----------------------------------------------------

    #calculate fermi energy Efi of intrinsic sc
    #WORKS
    Nc_300k = 2.78e+25  # 1/m3
    Nv_300k = 9.84e+24  # 1/m3
    Eg_300k = 1.12 * _ev_to_J  # J
    T2 = 287

    efi = fermi_ernergy__intr_C(Nc_300k,Nv_300k,Eg_300k,T2)
    print("EFi_SC = ",efi/_ev_to_J," ev")

    print(2*_kb*300/fermi_energy_metal_0K(ne))






main()