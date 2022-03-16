# Muehleder Christoph 01604413
import numpy as np
import matplotlib.pyplot as plt

# read file
data_pre = np.loadtxt("reaction_rate.txt",delimiter=',')
data_pre2 = np.transpose(data_pre)

# create important variables
S = data_pre2[0]                # substrate concentration in percent
V = data_pre2[1]                # reaction velocity
N_data = len(data_pre2[0])      # number of data points
n_fit = 2                       # a1*1, a2*1/s
beta = 1/V                      # beta = 1/v1 , 1/v2 , 1/v3...

# print data
#print("S[] = ")
#print(S)
#print("\n\n")
#rint("V[] = ")
#print(V)

# create matrix and fill
A = [[0 for n in range(n_fit)] for m in range(N_data)]  # A is m x n .. m rows, n cols

for i in range(N_data):
    for j in range(n_fit):
        A[i][j] = np.power(S[i], -j)        # 1, 1/s1,
                                            # 1, 1/s2

A_inv = np.linalg.pinv(A)                   # pseudo inverted matrix
a_j = A_inv.dot(beta)                        # get fit parameters aj = A_inv * beta ; b
print(a_j)
# a[0] = a
# a[1] = b
# (S) = a + b/S
# stdev of fitting parameters.. 3.4.1 ??

# create new data for reaction velocity and its inverse function
f1_new = np.dot(A,a_j)                      # f = beta = 1/V_j = A*a
V1_new = 1/f1_new                           # V_jj = 1/beta_j

# plot reaction velocity
plt.figure()
plt.plot(S, V)
plt.plot(S, V1_new)
plt.ylabel("V")
plt.title("Reaction Velocity")
plt.legend(["File Data","Fit Function"])
plt.xlabel("[S] / %")
plt.show()

# plot inverse reaction velocity
plt.figure()
plt.plot(S, 1/V)
plt.plot(S, f1_new)
plt.ylabel("1/V ")
plt.title("Inverse Reaction Velocity")
plt.legend(["File Data","Fit Function"])
plt.xlabel("[S] / %")
plt.show()