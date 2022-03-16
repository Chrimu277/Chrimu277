import random

from datetime import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# 1) create sine wave emitter
# 2) let wave travel forth against target
# 3) build reflection
# 4) interfere

A0 = 1
phi0 = 0
t0 = 0
_lambd = 500e-9
_c = 3e+8


def A(A0:float,w:float,t:float,k:float,x:float,phi0:float):
    return A0*np.cos(k*x-w*t+phi0)

# wave, t=0, stationary, nothing moves
def test1():
    x_arr = np.arange(0, 1e-6, 1e-8)  # m
    E_arr = []
    for i in range(len(x_arr)):
        E_arr.append(A(A0,w=2*np.pi*_c/_lambd,t = 0, k= 2*np.pi/_lambd, x = x_arr[i], phi0 = 0))

    print(E_arr)
    plt.plot(x_arr,E_arr)
    plt.show()

# wave travels forth
def test2():

    timesteps = 1000 # 1000 times
    positions = 1000 #1000 steps

    dx = 1e-9 #10 nm
    dt = 3e-17 #1 us

    t_arr = np.arange(0,timesteps*dt,dt)
    x_arr = np.arange(0,positions*dx,dx)

    # save results
    E_matrix = []

    fig, ax = plt.subplots()

    for j in range(timesteps):
        print("t = ",t_arr[j])
        E_t = np.zeros((len(x_arr)))
        for i in range(len(x_arr)):
            E_t[i] = A(A0=1,w=2*np.pi*_c/_lambd,t=t_arr[j],k=2*np.pi/_lambd,x= x_arr[i],phi0=0)
        ax.clear()
        ax.plot(x_arr,E_t)
        plt.xlabel("x-direction [m]")
        plt.ylabel("E-Field strength [V/m]")
        plt.pause(0.000001)

#wave gets emitted, reflected, bounces back + interferes with itself

class E_wave_particle:

    E0 = 1
    x = np.zeros(3,dtype=float)
    dir = np.asarray([1,1,1],dtype=float)
    c = 1

    def __init__(self):
        pass

    def propagate(self,dx:list):
        self.x[0] += dx[0]*self.dir[0]
        self.x[1] += dx[1]*self.dir[1]
        self.x[2] += dx[2]*self.dir[2]

    def reflect(self,axis=0):
        if axis > 2:
            return
        self.dir[axis] = -self.dir[axis]



def test3():

    timesteps = 100 # 1000 times
    positions = 20 #1000 steps

    dx = 1
    dt = 1

    t_arr = np.arange(0,timesteps*dt,dt)
    x_arr = np.arange(0,positions*dx,dx)

    # save results
    photons = []

    fig, ax = plt.subplots()

    N_p = 2

    for j in range(int(len(t_arr))):
        E_TOT = np.zeros(len(x_arr))
        # local field strength

        if len(photons) < N_p:
            E_t =  A(A0=1,w=2*np.pi*_c/_lambd,t=t_arr[j]*1e-2,k=2*np.pi/_lambd,x= 0,phi0=0)
            particle = E_wave_particle()
            particle.E0 = E_t
            particle.x[0] = np.random.randint(0,positions)
            photons.append(particle)
            print("created particle at ",particle.x)


        for i in range(len(photons)):
            photons[i].propagate([dx,0,0])
            print(photons[i].x)
            if photons[i].x[0] >= x_arr[-1] or photons[i].x[0]<= x_arr[0]:
                photons[i].reflect(0)
        print("--")

        for p in photons:
            print("x-coordinate: ",p.x[0])
            E_TOT[int(p.x[0])] += p.E0

        print(E_TOT)

        ax.clear()
        ax.plot(x_arr,E_TOT)
        plt.xlabel("x-direction [m]")
        plt.ylabel("E-Field strength [V/m]")
        #plt.ylim((-1,1))
        plt.pause(0.000001)


if __name__ == '__main__':
    test3()