import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import tifffile

def In(t,A,B,C,D,t0):
    theta = (t+t0)*np.pi/180
    return 0.5*(A+ B*np.sin(2*theta) + C*np.cos(2*theta) + D*np.sin(4*theta))


def A(I, theta):
    if len(I) != len(theta):
        return -1
    A = 0.0
    for i in range(len(I)):
        A += I[i]
    return A * 2/len(I)
def B(I, theta):

    B = 0.0

    if len(I) != len(theta):
        return -1

    for i in range(len(I)):
        B += I[i]*np.sin(2*theta[i]*np.pi/180)

    return B * 4/len(I)
def C(I, theta):
    C = 0.0
    if len(I) != len(theta):
        return -1
    if len(I) != len(theta):
        return -1

    for i in range(len(I)):
        C += I[i] * np.cos(4 * theta[i]*np.pi/180)
    return C * 4/len(I)
def D(I, theta):

    D = 0.0

    if len(I) != len(theta):
        return -1

    for i in range(len(I)):
        D += I[i] * np.sin(4 * theta[i]*np.pi/180)

    return D * 4/len(I)

def stokesvektor(a:float,b:float,c:float,d:float,norm = False):
    sv = np.array([a - c, 2 * c, 2 * d, b],dtype=float)
    if norm:
        sv = sv / sv[0]
    return sv

def norm(vec):
    norm = 0.0
    for i in range(len(vec)):
        norm += vec[i]**2

    return np.sqrt(norm)


basename = "source/stokes_measure_lobes_deg"
index = [str(i) for i in np.arange(0,180,10)]
filetype = ".tiff"

# angles 0 .. 170 deg + theta0
theta0 = [8]
theta_n1 = np.arange(0,180,10)

intensity_matrices = [] # saves intensity of each image and pixel [i][l][h]
L = 0
H = 0


# GET IMAGES
for i in range(len(theta_n1)):      # for each image

    im = tifffile.imread(basename+index[i]+filetype) # load image

    intensity_matrix = np.asarray(im,dtype=float)   # cast to float matrix

    L = len(intensity_matrix)                       # read L x H
    H = len(intensity_matrix[0])

    maxI = np.max(intensity_matrix)                 # norm to 256 intensity
    for k in range(len(intensity_matrix)):
        for j in range(len(intensity_matrix[0])):
            intensity_matrix[k][j] *= 256/maxI

    plt.imshow(intensity_matrix, interpolation='nearest')
    #plt.show()
    plt.close()

    print("loading image with size", intensity_matrix.shape)
    intensity_matrices.append(intensity_matrix)


# select window
L0 = 0
H0 = 0
# 0, 0
dL = 768
dH = 1024
# 768, 1024

print("loaded all images")
stokes_all = np.zeros((dL,dH,5),dtype=float)
print("len matrices: ",len(intensity_matrices))


# GET STOKES VECTOR
totalsteps = dL*dH

for n in range(len(theta0)):

    theta_n = theta_n1 + theta0[n]
    print("calc with .. ",theta_n)

    donesteps = 0
    for l in range(dL):          # for every pixel
        for m in range(dH):
            donesteps += 1

            #print("l, m = ",l,m)
            Ilm = [intensity_matrices[i][L0+l][H0+m] for i in range(len(intensity_matrices))]

            A1 = A(Ilm, theta_n)
            B1 = B(Ilm, theta_n)
            C1 = C(Ilm, theta_n)
            D1 = D(Ilm, theta_n)
            s = stokesvektor(A1,B1,C1,D1, norm = False)

            stokes_all[l][m][0:4] = s
            stokes_all[l][m][4] = 0



    # make matrices of stokes vector
    S0_matrix = np.zeros((dL,dH))
    S1_matrix = np.zeros((dL,dH))
    S2_matrix = np.zeros((dL,dH))
    S3_matrix = np.zeros((dL,dH))
    P_matrix = np.zeros((dL, dH))

    for l in range(dL):
        for h in range(dH):
            S0_matrix[l][h] = stokes_all[l][h][0]
            S1_matrix[l][h] = stokes_all[l][h][1]
            S2_matrix[l][h] = stokes_all[l][h][2]
            S3_matrix[l][h] = stokes_all[l][h][3]


    # norm matrix
    maxS0 = np.max(S0_matrix)
    print("min max S0", np.min(S0_matrix), " ", np.max(S0_matrix))
    print("min max S1", np.min(S1_matrix), " ", np.max(S1_matrix))
    print("min max S2", np.min(S2_matrix), " ", np.max(S2_matrix))
    print("min max S3", np.min(S3_matrix), " ", np.max(S3_matrix))
    S0_matrix *= 1 / maxS0
    S1_matrix *= 1 / maxS0
    S2_matrix *= 1 / maxS0
    S3_matrix *= 1 / maxS0


    # save results as picture
    resultname = "results/"
    fig,axes = plt.subplots(nrows=2, ncols=2)
    plt.suptitle("Stokes Parameters normed to 1")

    im1 = axes[0][0].imshow(S0_matrix,vmin=-1,vmax=1,cmap='jet')#,vmin=-1,vmax=1)
    axes[0,0].set_title("S0")
    plt.colorbar(im1)

    im2 = axes[0,1].imshow(S1_matrix,vmin=-1,vmax=1,cmap='jet')  # ,vmin=-1,vmax=1)
    axes[0,1].set_title("S1")

    im3 = axes[1, 0].imshow(S2_matrix,vmin=-1,vmax=1,cmap='jet')  # ,vmin=-1,vmax=1)
    axes[1, 0].set_title("S2")

    im4 = axes[1, 1].imshow(S3_matrix,vmin=-1,vmax=1,cmap='jet')  # ,vmin=-1,vmax=1)
    axes[1, 1].set_title("S3")


    plt.savefig(resultname + f"Stokes_{theta0[n]}.jpg")
    plt.show()
    plt.close()
