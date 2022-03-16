# Muehleder Christoph 01604413
import numpy as np
import matplotlib.pyplot as plt
import time


# ----FUNCTIONS DEFINITIONS--------


# gets pixel value from (x,y) coordinate values
def getPixelValue_x_y(x, y):
    if np.abs(x) < (_COLCOUNT / 2 - 1) and np.abs(y) < (_ROWCOUNT / 2 - 1 ):
        arr_y = int(np.rint(_ROWCOUNT / 2 - y))
        arr_x = int(np.rint(_COLCOUNT / 2 + x))
        return _FILEDATA[arr_y][arr_x]
    else:
        return 0


# maps (xi,eta,phi) to (x,y)
def getPixelValue_xi_eta_phi(xi, eta, phi):
    x = np.cos(phi) * xi - np.sin(phi) * eta                            # Assignment formula
    y = np.sin(phi) * xi + np.cos(phi) * eta
    return getPixelValue_x_y(x, y)


# goes through image on fixed xi to create a projection value
def p_phi_xi(phi, xi):
    sum = 0.0
    eta = -_RADIUS
    while eta <= _RADIUS:
        sum += _STEPSIZE_ETA * getPixelValue_xi_eta_phi(xi, eta, phi)
        eta += _STEPSIZE_ETA
    return sum


# ... row col -> x y
def back_integration(col, row, projections):
    x = col - _COLCOUNT/2
    y = _ROWCOUNT/2 - row
    sum = 0.0

    for i in range(len(_PHIS_DEG)):
        phi_rad1 = _PHIS_DEG[i] * np.pi / 180                       # degrees to rad
        xi_dash = get_XI(x, y, phi_rad1)                            # get xi' from x,y,phi
        xi_position = int(len(_XIS) / 2 + (xi_dash/_STEPSIZE_XI))   # gets position of xi' in _XIS[]
        sum += projections[i][xi_position]  # p~ [phi][xi]          # sums up projections[phi][xi]

    return sum / (2*(len(_PHIS_DEG)+1))


# returns xi from x,y,phi
def get_XI(x, y, phi):
    xi : float = x * np.cos(phi) + y * np.sin(phi)                  # //Assignment formula
    return xi


# ------PROGRAM START------------------------------------------------

# Read File data---

_FILEDATA = []
file = open("bird.txt", 'r')                                # open file

for line in file:                                           # for each line
    tokens = line.split(',')                              # split pixel values by comma
    _FILEDATA.append(np.asarray(tokens).astype(np.int))   # converts tokens-list to tokens_array, converts to int array

file.close()                                                # close file

# show bird
plt.figure()
plt.title("Original File Data")
plt.imshow(_FILEDATA, cmap=plt.cm.binary)
plt.show()

# create all important variables

_COLCOUNT = len(_FILEDATA[0])                                               # nr of columns
_ROWCOUNT = len(_FILEDATA)                                                  # nr of rows
_RADIUS = int(np.round(np.sqrt(_COLCOUNT ** 2 + _ROWCOUNT ** 2) * 1 / 2))   # maximum xi value (radius of circle)
_STEPSIZE_ETA = 3                                                           # delta eta
_STEPSIZE_XI = 3                                                            # delta xi
_STEPSIZE_PHI = 3                                                           # delta phi
_PHIS_DEG = list(range(0, 180, _STEPSIZE_PHI))                              # contains every angle in degrees
_XIS = list(range(-_RADIUS, _RADIUS, _STEPSIZE_XI))                         # contains every xi value

# create k_xi[j]
k_xij = np.zeros(len(_XIS))                                 # k_xi,j
k_xij_opt = np.fft.fftfreq(len(_XIS))                       # returns optimized kxij values
for i in range(len(_XIS)):
    k_xij[i] = (2*np.pi* i) / len(_XIS)*_STEPSIZE_XI        # k_xi_j = 2pi j / (Nxi * delta)    //Assignment formula

# sweep through image through all angles
# Make projections through image
projections_data = [[0 for i in range(len(_XIS))] for j in range(len(_PHIS_DEG))]

start1 = time.time()                                                # measure time

for i in range(len(_PHIS_DEG)):                                     # for every angle
    phi_rad = _PHIS_DEG[i] * (np.pi / 180);
    for j in range(len(_XIS)):                                      # for every xi
        projections_data[i][j] = p_phi_xi(phi_rad, _XIS[j])         # create projection of (phi,xi)
        print(str(_PHIS_DEG[i]) + "°/ XI = " + str(_XIS[j]))

end1 = time.time()                                                  # measure time


# Surface Plot of projections

# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# X1= np.asarray(_XIS)
# Y1 = np.asarray(_PHIS_DEG)
# X1, Y1 = np.meshgrid(X1, Y1)
# Z1 = np.asarray(projections_data)
#
# plt.figure()
# surf = ax.plot_surface(X1, Y1, Z1.reshape(X1.shape))
# plt.xlabel("xi")
# plt.ylabel("phi / rad")
# plt.show()

# Heat Plot of Projections
plt.figure()
p = plt.imshow(projections_data)
plt.colorbar(p)
plt.show()

# --c---2d Fourier transformation of projections

projections_fourier = projections_data.copy()

for i in range(len(_PHIS_DEG)):  # for each phi, make fourier transformation (projection slice theorem)
    fou = np.fft.fft(projections_data[i]);
    for j in range(len(_XIS)):
        projections_fourier[i][j] = fou[j]

# power density plot of 1 projection slice at phi = 0

plt.figure()
plt.plot(k_xij,(np.abs(projections_fourier[0])/len(_XIS))**2)
plt.title("Power density Plot of Fourier Transformed Porjections phi = 0")
plt.show();

# ---e--- show 1 projection slice at fixed angle phi

plt.figure()
plt.plot(_XIS, projections_data[0])
plt.title("Projection at phi = "+str(_PHIS_DEG[0]))
plt.xlabel("Xi")
plt.ylabel("Sum of Density (projection)")
plt.show()

# ---g--- modified back projection, inverse fft
projections_fourier_inverse = [[0 for i in range(len(_XIS))] for j in range(len(_PHIS_DEG))]  # p~ (xi)
projections_fourier_optimized = projections_fourier.copy()
projections_fourier_inverse_optimized = projections_fourier_inverse.copy()

for i in range(len(_PHIS_DEG)):
    for j in range(len(_XIS)):
        projections_fourier[i][j] *= np.abs(k_xij[j])  # multiply k_xi on F(k_xi*cos phi, k_xi*sin phi)
        projections_fourier_optimized[i][j] *= np.abs(k_xij_opt[j])

    ifou = np.fft.ifft(projections_fourier[i]);
    ifou_opt = np.fft.ifft(projections_fourier_optimized[i]);

    for j in range(len(_XIS)):
        projections_fourier_inverse[i][j] = np.real(ifou[j])
        projections_fourier_inverse_optimized[i][j] = np.real(ifou_opt[j])

plt.figure()
plt.title("Inverse fourier transformed projections")
p = plt.imshow(projections_fourier_inverse)
plt.colorbar(p)
plt.show()


# ---h--- discretize real space, find xi' , get f(x',y')
generated_image_xy = np.zeros(shape=(_ROWCOUNT, _COLCOUNT))  # picture [y][x]
generated_image_xy_opt = np.zeros(shape=(_ROWCOUNT, _COLCOUNT))

start = time.time()
for y in range(_ROWCOUNT):
    print("generate line " + str(y))
    for x in range(_COLCOUNT):
        generated_image_xy[y][x] = back_integration(x,y,projections_fourier_inverse)
        generated_image_xy_opt[y][x] = back_integration(x, y, projections_fourier_inverse_optimized)

end = time.time()
print("duration for generating pic: "+ str(end-start)+" sec")

# show image without
p = plt.imshow(generated_image_xy, cmap=plt.cm.binary)
plt.colorbar(p)
plt.title("Normal")
plt.show()

p = plt.imshow(generated_image_xy_opt, cmap=plt.cm.binary)
plt.colorbar(p)
plt.title("Frequency optimized")
plt.show()
print("duration for generating projections: " + str(end1-start1)+" sec")