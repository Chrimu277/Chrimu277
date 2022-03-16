
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
from scipy.stats import randint
import math as mat
from datetime import datetime
from sympy import *

def stdev_N(Ni,N):
    return np.sqrt(Ni*(1-(Ni/N)))

def stdev_N_simple(Ni):
    return np.sqrt(Ni)

def stdev_p(Ni,N):
    return (1/N) * stdev_N(Ni,N)

def toPercent(p):
    return str(p*100)+" %"

# A) how good are these numbers?
print(" --- A) ---------------------------------")
N = 35.0
pi = 0.938 #p_i

# i) How many tests were FALSE NEGATIVE.. ANSWER: N_falneg
Ni = np.round(pi*N,0) #N_i
Nf = np.round(N - Ni,0) #N_f

print("of ",str(N)," positive people "+str(Ni)+" were tested true positive")
print("so ",str(Nf)," were tested false negative")

# ii) What is uncertainty on that number?
sig_Ni = stdev_N_simple(Ni)
sig_pi = sig_Ni/N
sig_Nf = stdev_N_simple(Nf)
sig_pf = sig_Nf/N

print("sni: +-",sig_Ni)
print("spi: +-",toPercent(sig_pi))

print("snf: +-",sig_Nf)
print("spf: +-",toPercent(sig_pf))

# How could one sensibly specify the actual result of study?
# iii) assume 94% as sensitivity. How many covid postiive
# people wold have had to be tested for STDEV of 0.1% or 1%?

pi = .94
sig_pi_w1 = 0.001 # sigma pi wanted 1 .. 0.1%
sig_pi_w2 = 0.01 #1%

N_test = 1

while True:
    sig_pi = np.sqrt(pi*N_test*(1-pi))/N_test
    if sig_pi < sig_pi_w2:
        break;
    N_test += 1

print(N_test," people have to be tested for sigma_ni < ",sig_pi_w2)

while True:
    sig_pi = np.sqrt(pi*N_test*(1-pi))/N_test
    if sig_pi < sig_pi_w1:
        break;
    N_test += 1

print(N_test," people have to be tested for sigma_ni < ",sig_pi_w1)

# B) Calculate ”running averages” in O(N) computer-time
# time series xi = 1,2,....


# a) write S(t) recursion in t
# should take O(N) steps

def S(t,xj):
    #print("call S(",t,")")
    N = len(xj)
    if(t<N-1):
        #print("= xj[",N-t,"] + S(",t+1,")\n")
        return xj[N-t-1]+S(t+1,xj)
    else:
        #print("= xj[0]\n")
        return xj[0]

# b) test
print(" --- b) ---------------------------------")
N = 10
xj = np.arange(N) # xj = [0,1,2..,8,9]
print("xj is ",xj)

arr_t = np.arange(N)
print("running b with t[] = ",arr_t)
for c in arr_t:
    print("----- t = ",c,"------")
    print("S(",c,") = ",S(c,xj))


# c) run program with N = 3 with xj = 1,2,3
# verify by hand
print("---- c) ---------")
xj = np.array([1,2,3])
arr_t = [0,1,2]

for t in arr_t:
    print("S(",t,") = ", S(t,xj))

# d) run program with several large values of N !
# compare time stamps
print("---- d) ---------")
Ns= np.linspace(10,750,100)
runtimes = []

for N in Ns: # i from 0 to 3

    # make array xj that is N long and filled with pi
    xj = []
    for j in range(0,int(N)):
        xj.append(np.pi)

    timestamp1 = datetime.now()

    #cycle through S(t) with each t inside (0,N)
    # should take O(N^2) then
    for t in range(0,int(N)):
        S(t,xj)

    timestamp2 = datetime.now()
    tdelta= timestamp2-timestamp1
    runtimes.append(tdelta.microseconds)

# plot runtimes
plt.figure()
plt.plot(Ns,runtimes)
plt.xlabel("N")
plt.ylabel("runtime in us")
plt.show()

