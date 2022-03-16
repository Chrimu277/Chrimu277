#Lineare Regression
import numpy as np
import matplotlib.pyplot as plt

#sqrd_mean: x^2 strich -----
def sqrd_mean(list):
    s = 0.0;
    for i in range(len(list)):
        s += np.power(list[i],2);

    s /= len(list);
    return s;

#fehlerniveau
s = 0.75

#load file
data = np.loadtxt('MCMC_UE2019.dat')

#x, y
x_list = data[0];
y_list = data[1];

#plot data
plt.scatter(x_list, y_list,25);
plt.title('MCMC_UE2019 Messdaten')
plt.xlabel('X');
plt.ylabel('Y');
plt.grid(True,'major','both');
plt.errorbar(x_list,y_list,s,0);
plt.show();

#Maximum likelihood : a,b

yavg = np.average(y_list);
xavg = np.average(x_list);

x2avg = sqrd_mean(x_list);
y2avg = sqrd_mean(y_list);

#auf blatt papier errechnet
a = yavg;
b = 0;

def f(a,x,b):
    return a*x+b;

#plot function

