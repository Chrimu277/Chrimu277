import numpy as np

a0 = 2.91       #angstrom

a1 = [1.0   ,1.0,   -1.0]
a2 = [-1.0  ,1.0,   1.0]
a3 = [1.0   ,-1.0,  1.0]


b1 = 2*np.pi*(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))
b2 = 2*np.pi*(np.cross(a3,a1))/(np.dot(a1,np.cross(a2,a3)))
b3 = 2*np.pi*(np.cross(a1,a2))/(np.dot(a1,np.cross(a2,a3)))

print("b1 is",b1)
print("b2 is",b2)
print("b3 is",b3)