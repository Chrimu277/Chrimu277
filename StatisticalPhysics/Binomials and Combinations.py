import random
from datetime import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spspe

# combinations
# calculate, how many combinations of k out of N there are
# order not important
a1 = spspe.comb(N=5,k=5, exact = False, repetition= False) # = binomial coefficient (N k)
print("Comb: 10 different; choose 5; no repetitons -> ",a1)
# is equal to binom(n,k)
# for N = k -> 1 option

a2 = spspe.comb(N=10,k=5, exact = False, repetition= True)
print("Comb: 10 different; choose 5; with repetitons -> ",a2)

# -------------

# permutations
# calculate how many permutations of j and k elements there are, N = j+k
# order matters
a3 = spspe.perm(N=4,k=2,exact=False)
print("Perm: 4 different; choose 2 at a time: -> ",a3)


# -----------------
# telephone numbers with N = 5, composed of 0 and 1s   n 0s  and m 1s
# for example: 01011; 01110
n = 2
m = 2
# total phone numbers possible (permutations)
n1 = spspe.perm(n+m,n+m)
print(n1)