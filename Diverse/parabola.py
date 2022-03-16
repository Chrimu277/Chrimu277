# Autor Mühleder Christoph
# April 4 2021
# DONE
# Schiefer Wurf / Parabel

import numpy as np
import matplotlib.pyplot as plt

dt = 0.001  #s
x0 = [0,0]  #m
g = -9.81   #m/s^2

v0 = [1,20]   # x-speed, y-speed
t0 = 0
t1 = 10

steps =  int((t1-t0) / dt)

x = [[],[]]     # speichert alle x,y Positionen
v = [[],[]]

x[0].append(x0[0])
x[1].append(x0[1])

v[0].append(v0[0])
v[1].append(v0[1])
t=[]
t.append(0)

for i in range(1,steps):
    t.append(dt*i)
    x[0].append(0)
    x[1].append(0)

    v[0].append(0)
    v[1].append(0)

    v[0][i] = v[0][i-1]
    v[1][i] = v[1][i-1] + g*dt

    x[0][i] = x[0][i-1] + v[0][i]*dt
    x[1][i] = x[1][i-1] + v[1][i]*dt


plt.plot(t,x[0])
plt.plot(t, x[1])
plt.show()

plt.plot(t,v[0])
plt.plot(t, v[1])
plt.show()