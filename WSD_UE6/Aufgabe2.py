#Spectral Analysis

#------------------------------------------------------- a) ----------------------------------------------------------

import math
import numpy as np
import matplotlib.pyplot as plt

N1 = 63;
N2 = 511;

zeitreihe1 = []
zeitreihe2 = []
zeitpunkte = []

A1 = 0.8
B1 = 0
f1 = 0.1
s1 = 1

A2 = 0
B2 = 10
f2 = 0.1
s2 = 1

for i in range(N1):
    #sekunden liste
    zeitpunkte.append(i);
    #random normalverteiltes rauschen erzeugen
    ei1 = np.random.normal(0, 1);
    ei2 = np.random.normal(0, 1);
    #zeitreihe erzeugen
    di1 = A1*math.cos(2*np.pi*f1*i)+B1*np.sin(2*np.pi*f1*i);
    di2 = A2*math.cos(2*np.pi*f1*i)+B2*np.sin(2*np.pi*f2*i);
    #werte addieren und zur zeitreihe
    zeitreihe1.append(di1+ei1);
    zeitreihe2.append(di2+ei2);

#plot signal 1 mit errorbar sigma1 s1
plt.plot(zeitpunkte,zeitreihe1,'k');
plt.title("Signal 1")
plt.errorbar(zeitpunkte,zeitreihe1,s1,0,'','r');
plt.xlabel("t [s]")
plt.ylabel("d(t)")
plt.show();

#plot signal 2 mit errorbar sigma2 s2
plt.plot(zeitpunkte,zeitreihe2,'k');
plt.title("Signal 2")
plt.errorbar(zeitpunkte,zeitreihe2,s2,0,'','r');
plt.xlabel("t [s]")
plt.ylabel("d(t)")
plt.show();

#--------------------------------------------------- b) ---------------------------------------------------------------

#zeitpunkte bis 511 mit sekunden auffüllen, schwingungen mit 0
for i in range(N1,N2):
    zeitpunkte.append(i);
    zeitreihe1.append(0);
    zeitreihe2.append(0);

#complex arrays mit fast fourier transformation (signal in spektrum zerlegen)
fft1 = np.fft.fft(zeitreihe1); #len = 511
fft2 = np.fft.fft(zeitreihe2);

#leistungsdichte |Sf)|    (betrag der fourierreihe : 1/N2 * (summe( re^2 + imag^2 )
leist_dichte1 = [];
leist_dichte2 = [];

#leistungsdichte-array aus fft-array errechnen und zeitpunkte und arrays auf 511s auffüllen
for i in range(N2):
        leist_dichte1.append((np.power(fft1[i].real, 2) + np.power(fft1[i].imag, 2))/N2);
        leist_dichte2.append((np.power(fft2[i].real, 2) + np.power(fft2[i].imag, 2))/N2);

#plot leistungsdichten
plt.plot(zeitpunkte, leist_dichte1);
plt.xlabel('f');
plt.ylabel('S(f)')
plt.title('FFT Signal 1')
plt.show();
#(man sieht, das rauschen führt durch die niedrige schwingungsamplitude bei f1 zu unklarem spektrum)
#bei S2(f2) ist die amplitude ca. 100* größer als rauschen, daher 2 klare peaks !
plt.plot(zeitpunkte, leist_dichte2);
plt.title('FFT Signal 2')
plt.xlabel('f');
plt.ylabel('S(f)')
plt.show();

#--------------------------------------------------- c) -------------------------------

d1_square_mean = 0;
d2_square_mean = 0;

#p(f|d,t,B)
p1 = []
p2 = []

for i in range(len(zeitreihe1)):
    d1_square_mean += zeitreihe1[i]**2;
for i in range(len(zeitreihe2)):
    d2_square_mean += zeitreihe2[i]**2;

for i in range(len(leist_dichte1)):
    p1.append(np.power(1 - ((2*leist_dichte1[i])/(N2*d1_square_mean)),(1-N2/2))) # formel aus angabe
for i in range(len(leist_dichte2)):
    p2.append(np.power(1 - ((2 * leist_dichte2[i]) / (N2 * d2_square_mean)), (1 - N2 / 2)))

plt.plot(zeitpunkte, p1 );
plt.xlabel('f');
plt.ylabel('p1(f|d,t,B)')
plt.show()

plt.plot(zeitpunkte, p2);
plt.xlabel('f');
plt.ylabel('p2(f|d,t,B)')
plt.show()

#p(f|d,t,B) proportional zur spektralanalyse
