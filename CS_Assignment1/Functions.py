# AUTHOR: Christoph Muehleder 01604413
# DATE: March 15, 2021
# Computer Simulations, Assignment 1, HELPER CLASS

import matplotlib.pyplot as plt
import numpy as np


class MyFunctionHelperClass:

    def __init__(self):
        pass

    # the Probability density function of exercise 3
    @staticmethod
    def pdf_g(x=None, z=0.976):
        if x is None:
            return 0;
        t1 = 1 / z;
        t2 = np.power(np.sin(x), 2)
        t3 = np.exp(-x)
        t4 = np.power((1 + np.exp(-x)), 2)
        return t1 * t2 * 2 * t3 / t4

    # the logistic distribution shall be envelope of g(x) (Exercise 3)
    @staticmethod
    def envelope_func1(x, arr=np.array([1, 1, -5, 0.01])):
        c = arr[0]
        k = arr[1]
        x0 = arr[2]
        offs = arr[3]

        t1 = (1 + np.exp(-k * (x - x0)))
        return c / t1 + offs;

    # the new envelope with 1/x2 behaviour (Exercise 3)
    @staticmethod
    def envelope_func2(x, arr=np.array([0, 1])):
        x0 = arr[0]
        z = arr[1]
        if isinstance(x,np.ndarray):
            x_helper = []
            for xi in x:
                if xi-x0 < 0:
                    x_helper.append(z);
                else:
                    trial = z/np.power((xi-x0),2)
                    if trial > z:
                        x_helper.append(z)
                    else:
                        x_helper.append(trial)
            return np.array(x_helper)
        else:
            if x < x0 or np.abs(x-x0) < 0.01:
                return z;
            else:
                trial = z / np.power((x - x0), 2)
                if trial > z:
                    return z;
                else:
                    return trial;

    # the logstic distribution of Exercise 2,  x ; element [0, inf] , h(x) = h(-x)
    @staticmethod
    def logistic_distr(x = 0.0):
        t1 = np.exp(-x)
        t2 = np.power((1+np.exp(-x)),2)
        return 2*t1/t2;

    @staticmethod
    def askUserInt(str,default = -1,min = 0,max = 0):
        ans = 0
        if min == 0 and max == 0:
            in_string = input(str)
            try:
                ans = int(in_string)
                return ans
            except:  # entered bullshit
                return default
        else:
            print("Value in [", min, " .. ", max, "]")
            in_string = input("")
            try:
                ans = int(in_string)
                return ans if (ans < max and ans > min) else default # check limits
            except:  # entered bullshit
                return default



