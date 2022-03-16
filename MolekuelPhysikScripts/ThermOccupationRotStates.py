import numpy as np
# not my source code

cm_in_J = 1.986e-23
eV_in_cm = 8065.5

w_e = 2358.57 # cm-1
w_ex_e = 14.34 # cm-1
B_e = 1.99824 # cm-1
alpha_e = 0.017318 # cm-1
D_e = 0.00000576 # cm-1

U_0 = 1.76 * eV_in_cm # eV
r_e = 1.09769 # A
T = 300.0 # K

h = 4.136e-15 * eV_in_cm # cm-1 s
c = 3e10 # cm/s
k_B = 0.695 # cm-1 / K

# rotational energy
def Erot (v,J):
    return h*c*(B_e - alpha_e*(v+0.5))*J*(J+1) + D_e * ((J*(J+1))**2)

def Evib (v):
    return h*c*w_e*(v+0.5) - h*c*w_ex_e*((v+0.5)**2)

# ratio of occupation probability
def P (i,E):
    return (2*i+1)*np.exp(-E/(k_B*T))

# Ground state
energy_ground = Erot(0,0) + Evib(0)
prob_ground = P(0,energy_ground)

print('Ground state:    Energy =',energy_ground,' Probability =',prob_ground)

# first 40 states
for i in range(40):
    energy = Erot(0,i) + Evib(0)
    probability = P(i,energy)
    print('State J =',i,':   E_J =',round(energy,7),'   P_J =',probability,'  P_J/P0 =',round(probability/prob_ground,3))