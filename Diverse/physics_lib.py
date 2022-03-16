
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
from scipy.stats import randint
import math as mat
from datetime import datetime

from sympy import *

#
# #define functions
# f_x = lambda x: x*x;
#
# def f_of_x(x):
#     return x*x;
#
# def my_sin(x):
#     return mat.sin(x);
#
# #simple integrate
# def integrate_array(arr,dx):
#     N = arr.size
#     S = np.sum(arr)*dx
#     return S;
#
# #simple integrate
# def integrate_function(fx,dx,x1,x2):
#     S = 0.5*(fx(x1)+fx(x2))
#     x= x1+dx
#     while x < x2:
#         S+= fx(x)
#         x+=dx
#
#     return S*dx
#
#
# # define x
# list1 = [1,2,3]
# N = int(10e6)
# dx = 0.001
# array1 = np.arange(0,2*np.pi,dx)
# linsp1 = np.linspace(0,100,N)
#
# # make f(x)
# y = f_x(array1)
# y2 = f_of_x(array1)
# siny = np.sin(array1)
#
# # show function
# plt.figure()
# plt.plot(array1,y)
# plt.show()
#
# #-----------------
#
# np.random.seed(datetime.now().microsecond)
# N = 2000
# i = np.random.rand(N) # gives N random numbers between 0 and 1
# j = i.copy()
# for index in range(len(j)):
#     j[len(j)-index-1] = j[len(j)-index-2]
#
# j[0] = 0
# print(i[0:5])
# print(j[0:5])
#
# correllation of xi and xi+1
# plt.scatter(i,j,s=.1)
# plt.show()
#
# # -------------------------
#
# j = np.arange(0,len(i))
# hist1,bins1 = np.histogram(i,bins=np.arange(0,1,0.1))
#
# print("histogram1: ",hist1)
# print("bins1: ",bins1)
#
# plt.figure()
# plt.hist(hist1,bins=bins1,range=[0,500],cumulative=False)
# plt.show()

np.random.seed(datetime.now().microsecond)

rand_2darr = np.random.random((3,2))
print(rand_2darr)