import numpy as np


#readfile
data_pre = np.loadtxt("reaction_rate.txt",delimiter=',')
data_pre2 = np.transpose(data_pre)

S = data_pre2[0] # substrate concentration in percent
V = data_pre2[1] # reaction velocity
N = len(data_pre2[0])

#V(S) = V0 / (K + S)
# f(S) = 1/V(S), fit parameters a   ... = K/V0  + S / V0
