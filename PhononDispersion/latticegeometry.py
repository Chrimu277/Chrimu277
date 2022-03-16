# ---- import globals ----
import numpy as np
# ---- import locals ----

# ---- FUNCTIONS ----

# returns 3 primitive lattice vectors of FCC
def primitiveLatticeVectorsFCC(a=1):
    a1 = 0.5 * a * (_x + _z)
    a2 = 0.5 * a * (_x + _y)
    a3 = 0.5 * a * (_y + _z)
    return a1,a2,a3

# returns 3 primitive lattice vectors of BCC
def primitiveLatticeVectorsBCC(a=1):
    a1 = 0.5 * a * (_x + _y - _z)
    a2 = 0.5 * a * (-_x + _y + _z)
    a3 = 0.5 * a * (_x - _y + _z)
    return a1,a2,a3

# squared absolute length of vector
def absVector(v:np.ndarray) -> float:
    return np.sqrt(np.dot(v,v))

# norms a vector to length 1
def normVector(v:np.ndarray) -> np.ndarray:
    return v / absVector(v)

# gives back 3 unit vectors (reciprocal) to a set 3 of input vectors
def reciprocalVectors(a1:np.ndarray,a2:np.ndarray,a3:np.ndarray):
    b1 = 2*np.pi * np.cross(a2,a3) / (np.dot(a1,np.cross(a2,a3)))
    b2 = 2 * np.pi * np.cross(a3, a1) / (np.dot(a2, np.cross(a3, a1)))
    b3 = 2 * np.pi * np.cross(a1, a2) / (np.dot(a3, np.cross(a1, a2)))
    return b1,b2,b3

# --- GLOBAL VARIABLES ----
_x = np.asarray([1, 0, 0])  # unit vectors in x,y,z
_y = np.asarray([0, 1, 0])
_z = np.asarray([0, 0, 1])

# ---- MAIN ----
if __name__ == '__main__':

    a = 1                       # lattice constant fcc
    a1,a2,a3 = primitiveLatticeVectorsFCC(a)
    print("FCC primitive lattice vectors:\n",np.matrix([a1,a2,a3]))

    b1, b2, b3 = primitiveLatticeVectorsBCC(a)
    print("BCC primitive lattice vectors:\n", np.matrix([b1,b2,b3]))

    c1,c2,c3 = reciprocalVectors(a1,a2,a3) # reciprocals of fcc
    print("rec. of fcc:\n",np.matrix([c1,c2,c3]))

    d1, d2, d3 = reciprocalVectors(b1, b2, b3) # rec. of bcc
    print("rec. of bcc:\n", np.matrix([d1,d2,d3]))