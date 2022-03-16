import numpy as np
import matplotlib.pyplot as plt

S1 = [1, -0.91, 0.01, -0.03]
S2 = [1,0,0.31,-0.76,0.811]

def In(t,A,B,C,D,t0):
    theta = (t+t0)*np.pi/180
    return 0.5*(A+B*np.sin(2*theta) + C*np.cos(2*theta) + D*np.sin(4*theta))

i1 = In()