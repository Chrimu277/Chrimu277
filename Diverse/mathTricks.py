
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
from scipy.stats import randint
import math as mat
from datetime import datetime

from sympy import *

def GoldenRatio():
    pass

def doMandelbrot():

    ct = int(3e3)
    maxit = int(1e2)

    y_min = -2
    y_max = 2

    x_min = -2
    x_max = 2

    xs = np.linspace(x_min,x_max,ct)
    ys = np.linspace(y_min, y_max,ct)
    grid = np.ndarray((len(xs),len(ys)),dtype=int)

    cnter = 0
    max = ct**2
    for i in range(ct):
        for j in range(ct):
            cnter = cnter+1
            c = complex(xs[i],ys[j])
            grid[i][j] = MandelbrotSteps(c,maxit)

    print(cnter)


    plt.figure(dpi=1000)
    plt.imshow(grid, cmap='Spectral', interpolation='bilinear')
    plt.xlabel("Re")
    plt.ylabel("Im")
    plt.show()

def MandelbrotSteps(C,max_iterations):
    z = C
    i = 0
    while np.abs(z) < 2 and i < max_iterations:
        z = z*z + C
        i = i + 1
    return i

def SquaredComplexNumbers():
    n = 10
    res = []
    ims = []
    np.random.seed(datetime.now().microsecond)
    c_start = complex(0.6, -0.6)
    leg = []

    for j in range(100):
        res.append([])
        ims.append([])
        im, re = np.random.rand(2)

        c_start = complex(re, im)

        for i in range(n):
            res[j].append(c_start.real)
            ims[j].append(c_start.imag)
            c_start = c_start * c_start
            if np.abs(c_start) > 2:
                break;

    fig = plt.figure()

    for j in range(10):
        plt.plot(res[j], ims[j])
        plt.ylim(-2, 2)
        plt.xlim(-2, 2)
        plt.show()

if __name__ == '__main__':

    # doMandelbrot()
    SquaredComplexNumbers()
    pass





