import numpy as np
import matplotlib.pyplot as plt

#function to integrate

#KONSTANTEN
h_bar = 1.054e-34
h = h_bar*2*np.math.pi
kb = 1.38e-23
n = 8.467e+28
w_D = 4.49e+13
w0 = 0
N = 1000
i = 0;
T = 200;
m_e = 9.109e-31

#stepsize
h = (w_D-w0)/N

#cv (T)
cv_ph_T = [];
cv_el_T = [];

#omegas (logarithmisch) und Ts aufsetzen
ws = np.arange(start = 1, stop = np.log(w_D), step = np.log(w_D)/N)
Ts = np.arange(start = 1, stop = 1e4, step = 10)

print(ws)

def func(w,T):
    fwt = np.power(w, 4) * np.exp(h_bar*w/(kb*T)) / (np.exp(h_bar*w/(kb*T)) - 1)**2
    return fwt

for j in range(len(Ts)):

    #cv_ph berechnen
    cv_ph_T.append(0.0)
    #numerisch integrieren
    for i in range(len(ws)-1):
        w1 = np.exp(ws[i])
        w2 = np.exp(ws[i+1])
        T1 = Ts[j]

        h = w2-w1

        fw1 = func(w1,T1)
        fw2 = func(w2,T1)

        cv_ph_T[j] += h * (fw1+fw2)/2

    cv_ph_T[j] *= 9*n*np.power(h_bar,5)/(np.power(kb,4)*np.power(343,3)*np.power(T1,2))

    #cv_ehl berechnen
    cv_el_T.append(0.0)
    cv_el_T[j] = np.power(np.pi/3,2/3)*m_e*np.power(n,1/3)*np.power(kb/h_bar,2)*T1




fig = plt.figure();
plt.plot(Ts,cv_ph_T)
plt.plot(Ts,cv_el_T)
plt.xscale('log')
plt.show()

#function for cv,el
# function for cv,ph
# plot -> ablesen