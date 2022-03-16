from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

def add(a,b):
    sum = a+b
    return sum

def mySin(x):
    ssq = np.sin(x)
    return ssq

def derivative(fx,x0 = 0,h=0.01):
    #derivate
    fx_x0 = (fx(x0+h)-fx(x0-h))/(2*h)
    return fx_x0

def myFunc(x):
    return x*x+5;

def rollStats():
    N_rolls = 4
    N_rep = 10
    dice_min = 1
    dice_max = 6

    sums = []

    for rep in range(N_rep):
        rolls = np.zeros(N_rolls)
        for i in range(N_rolls):
            rolls[i] = (np.random.randint(dice_min, dice_max, dtype=int))

        rolls_orderd = np.sort(rolls)
        highest3 = rolls_orderd[1:4]
        sumh3 = np.sum(highest3)
        sums.append(sumh3)
        if sumh3 > 10:
            print("good roll ! ",sumh3)
        else:
            print("oops ", sumh3)

    return sums

def rollRickey():
    N_rolls = 10
    dice_min = 1
    dice_max = 10
    limit = 8
    success = 0

    rolls = []
    for i in range(N_rolls):
        roll = (np.random.randint(dice_min, dice_max, dtype=int))
        rolls.append(roll)
        if roll >= limit:
            success = success + 1

    return success, rolls


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    succ, rolls = rollRickey()
    print("rolled ",succ, " successes")
    print("rolls:", rolls)