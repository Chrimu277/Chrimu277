import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import string
import array
from mpl_toolkits.mplot3d import Axes3D


# ----FUNCTIONS DEFINITIONS--------

# gets pixel value from (x,y) values
def func1(x, y):
    if np.abs(x) < _COLCOUNT / 2 - 1 and np.abs(y) < _ROWCOUNT / 2 - 1:
        arr_y = int(_ROWCOUNT / 2 + y)
        arr_x = int(_COLCOUNT / 2 + x)
        return _FILEDATA[arr_y][arr_x]
    else:
        return 0


# maps (xi,eta,phi) to (x,y)
def func2(xi, eta, phi):
    x = np.cos(phi) * xi - np.sin(phi) * eta
    y = np.cos(phi) * xi + np.cos(phi) * eta
    return func1(x, y)


# goes through image on fixed xi to create a projection value
def proj_phi(phi, xi, L, step):
    sum = 0.0
    eta = -L
    while eta < L:
        sum += step * func2(xi, eta, phi)
        eta += step
    return sum


def back_integration(col, row):
    # ... row col -> x y
    x = col - _COLCOUNT/2
    y = _ROWCOUNT/2 - row
    sum = 0.0
    for i in range(len(_PHIS_DEG)):
        phi_rad1 = _PHIS_DEG[i] * np.pi / 180;
        xi_dash = get_XI_dash(x, y, phi_rad1)
        xi_position = int(len(_XIS) / 2 + xi_dash/_STEPSIZE_XI)
        sum += projections_fourier_inverse[i][xi_position]  # p~ [phi][xi]
    return sum / (2*(len(_PHIS_DEG)+1))


def get_XI_dash(xdash, ydash, phi):
    xi_dash = xdash * np.cos(phi) + ydash * np.sin(phi)
    return np.round(xi_dash)


# ------PROGRAM START------------------------------------------------
#
# #Read in File data---
_FILEDATA = []
file = open("bird.txt", 'r')
for line in file:
    line_arr = line.split(',')
    val_arr = []
    for token in line_arr:
        val_arr.append(int(token))
    _FILEDATA.append(val_arr)

# show bird
#plt.imshow(_FILEDATA)
#plt.show()

# create important variables
_COLCOUNT = len(_FILEDATA[0])
_ROWCOUNT = len(_FILEDATA)
_RADIUS = int(np.round(np.sqrt(_COLCOUNT ** 2 + _ROWCOUNT ** 2) * 1 / 2)+1)  # L
_STEPSIZE_ETA = 2
_STEPSIZE_XI = 2
_STEPSIZE_PHI = 2

# Make projections through image
_PHIS_DEG = list(range(0, 180, _STEPSIZE_PHI))
_XIS = list(range(-_RADIUS, _RADIUS, _STEPSIZE_XI))
projection_data = [[0 for i in range(len(_XIS))] for j in range(len(_PHIS_DEG))]

print("xis length = " + str(len(_XIS)))
print(_XIS)
print("k_xi_j length = " + str(len(_PHIS_DEG)))
print(_PHIS_DEG)

# sweep through image thorugh all angles
for i in range(len(_PHIS_DEG)):
    phi_rad = _PHIS_DEG[i] * (np.pi / 180);
    for j in range(len(_XIS)):
        projection_data[i][j] = proj_phi(phi_rad, _XIS[j], _RADIUS, _STEPSIZE_ETA)
        print(str(_PHIS_DEG[i]) + "°/ XI = " + str(_XIS[j]))

# Surface Plot
fig = plt.figure()
ax = fig.gca(projection='3d')

X1= np.asarray(_XIS)
Y1 = np.asarray(_PHIS_DEG)
X1, Y1 = np.meshgrid(X1, Y1)
Z1 = np.asarray(projection_data)

surf = ax.plot_surface(X1, Y1, Z1.reshape(X1.shape))
plt.xlabel("xi")
plt.ylabel("phi / rad")
plt.show()

plt.imshow(Z1)
plt.xlabel("xi")
plt.ylabel("phi / rad")
plt.show()

# --c---2d Fourier transformation of projections
projections_fourier = [[0 for i in range(len(_XIS))] for j in range(len(_PHIS_DEG))]  # F[p(xi)] = F(k_xi..)
for i in range(len(_PHIS_DEG)):  # for each phi, make fourier transformation (projection slice theorem)
    fou = np.fft.fft(projection_data[i]);
    for j in range(len(_XIS)):
        projections_fourier[i][j] = fou[j]

print("proj fourier")
print(projections_fourier)
plt.figure()
p = plt.imshow(np.real(projections_fourier))
#plt.colorbar(p)
plt.show()

# create kxis
k_xij = np.zeros(len(_XIS))  # k_xi,j
print("k_xi_j length = " + str(len(k_xij)))
for i in range(len(_XIS)):
    k_xij[i] = (2*np.pi* i) / _STEPSIZE_XI
print("k_xij :" + str(k_xij))

# ---e--- show 1 projection slice at fixed angle phi

plt.figure()
plt.title("Fourier Transformed projection slice, phi = " + str(_PHIS_DEG[0]))
plt.xlabel("k_xi")
plt.ylabel("intensity")
plt.plot(k_xij, projections_fourier[0])
plt.show()

# ---g--- modified backprojection, inverse fft
projections_fourier_inverse = [[0 for i in range(len(_XIS))] for j in range(len(_PHIS_DEG))]  # p~ (xi)
for i in range(len(_PHIS_DEG)):
    for j in range(len(_XIS)):
        projections_fourier[i][j] *= np.abs(k_xij[j])  # multiply k_xi on F(k_xi*cos phi, k_xi*sin phi)

    ifou = np.fft.ifft(projections_fourier[i]);
    for j in range(len(_XIS)):
        projections_fourier_inverse[i][j] = np.real(ifou[j])

print("proj fourier inverse")
print(projections_fourier_inverse)
plt.figure()
p = plt.imshow(projections_fourier_inverse)
plt.colorbar(p)
plt.show()


# ---h--- discretize real space, find xi' , get f(x',y')
generated_image_xy = np.zeros(shape=(_ROWCOUNT, _COLCOUNT))  # picture [y][x]
for y in range(_ROWCOUNT):
    for x in range(_COLCOUNT):
        generated_image_xy[y][x] = back_integration(x,y)
        print("get "+str(y)+"/"+str(x))

p = plt.imshow(generated_image_xy)
plt.colorbar(p)
plt.show()
