import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#Pareto Verteilung

#-------------------------------------------------- a) ----------------------------------------------------------------

#warhscheinlichkeit paretoverteilung
def p(x,xmin,k):
    pi = 1;
    if(x < xmin):
        pi = 0;
    else:
        pi = np.power((xmin/x),(k+1)) * k/xmin;

    return pi;

#likelihood von messdaten xdata[]
def L(xdata,xmin,k):
    likeli = 1;
    for i in range(len(xdata)):
        likeli *= p(xdata[i],xmin,k);

    return likeli;


# generiere 20 Datenpunkte (xdata)
k = 2
x_min = 1
length = 20

x_data = [] #list

for i in range(20):
    x_data.append(np.random.pareto(k,x_min)[0])

#2x 20er list (achsen für plot)
ks = np.linspace(0,10,20); #x
x_mins = np.linspace(0.01, min(x_data), 20); #y

#20x20 matrix für p(k,xmins)... z
matrix = [[0 for x in range(length)] for y in range(length)] #type: list

#fülle matrix mit Likelihood funktionen auf
for i in range(length):
    for j in range(length):
        matrix[i][j] = L(x_data,x_mins[i],ks[j]);


#aus den listen ks und xmins wird --> meshgrid 2d für den plot
Xmesh,Ymesh = np.meshgrid(x_mins,ks);

#plot 3d
fig = plt.figure()
ax = Axes3D(fig);
ax.plot_surface(Xmesh,Ymesh,np.array(matrix));
ax.set_xlabel('x_min')
ax.set_ylabel('k')
ax.set_zlabel('P(x_data)')
plt.show()

#---------------------------------------------------- c) --------------------------------------------------------------

#erwartungswert und standardabweichung für k,xmin aus p(k,xmin|xdata,B)
