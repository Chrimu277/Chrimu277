
# AUTHOR: Christoph Muehleder 01604413
# DATE: March 13, 2021
# Computer Simulations, Assignment 1, Exercise 3

# --- global imports ---
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate as scipy_int
from datetime import datetime

# --- local imports ---
from CS_Assignment1.CS_A1_2 import InverseTransformation as Itm
from CS_Assignment1.Functions import MyFunctionHelperClass as MyFunctions

class RejectionMethod:

    # --- member variables ---

    # --- constructor ---
    def __init__(self):
        pass

    # rejection sampling of 2 pdf's
    @staticmethod
    def getSamplesRejectionMethod(N_samples, func_g, func_h,
                                  c, A, B, x_ct,
                                  args_h = None, args_g = None,
                                  plotit = False, printit = False, bunchsize = int(1e3)):

        # seed and make x_array
        np.random.seed(datetime.now().microsecond)
        x_arr = np.linspace(A, B, x_ct)
        dx = float((B - A) / x_ct)

        # PDF that we should follow & CDF
        g_x = func_g(x_arr) if (args_g is None) else func_g(x_arr,args_g)
        cdf_g_x = scipy_int.cumtrapz(g_x,x_arr,dx)
        norm_g = float(np.max(cdf_g_x))
        g_x = g_x / norm_g  # normed g(x) OK
        cdf_g_x = cdf_g_x / norm_g  # normed G(x)  OK

        # h(x) and H(x)
        h_x = func_h(x_arr) if (args_h is None) else func_h(x_arr,args_h)

        cdf_h_x = scipy_int.cumtrapz(h_x,x_arr,dx)
        norm_h = float(np.max(cdf_h_x))
        h_x = h_x / norm_h  # normed h(x)  OK
        cdf_h_x = cdf_h_x / norm_h  # normed h(x)  OK

        # plot envelope and pdf -- OK
        if printit: print(" ... plot g,G,h,h*c,H")
        if plotit:
            plt.figure()
            plt.plot(x_arr, g_x)
            plt.plot(x_arr[1:],cdf_g_x)
            plt.plot(x_arr, c * h_x)
            plt.plot(x_arr, h_x)
            plt.plot(x_arr[1:], cdf_h_x)
            plt.legend(["g(x)", "G(x)", "c*h(x) envelope", "h(x)", "H(x)"])
            plt.title("g, G, c*h, h, H")
            plt.xlabel("x")
            plt.ylabel("f(x)")
            plt.grid()
            plt.show()

        # generate x_T from h(x) through ITM sampling

        acc = []  # accepted values
        x_Ts = []  # sampled values from h(x)
        ctr = 0  # attempt counter

        if printit: print("drawing ",N_samples," samples from rejection method")

        while len(acc) < N_samples:

            x_T = Itm.getItmSampleFromPdf(func_pdf=func_h, A=A,B=B,x_ct=x_ct,
                N_sample=bunchsize, func_args = None if (args_h is None) else args_h)  # samples from h(x)
            r = np.random.rand(bunchsize)

            for i in range(bunchsize):
                h_xT = func_h(x_T[i]) if (args_h is None) else func_h(x_T[i], args_h)
                g_xT = func_g(x_T[i]) if (args_g is None) else func_g(x_T[i], args_g)

                h_xT = h_xT / norm_h
                g_xT = g_xT / norm_g

                ctr = ctr + 1
                x_Ts.append(x_T[i])
                if g_xT > r[i] * c * h_xT:
                    acc.append(x_T[i])


        # plot the sampled numbers
        if printit: print("plotting histogram of xT..")
        if plotit:
            plt.hist(x_Ts, bins=50, density=True)
            plt.title("Histogram of sampled x_T's of h(x)")
            plt.xlabel("x_T")
            plt.ylabel("density")
            plt.grid()
            plt.show()

        # acceptance rate and c
        acc_rate = (len(acc) / ctr)

        return acc[:N_samples],acc_rate


if __name__ == '__main__':

    # START
    print("\n--- executing 3) ---")
    t1 = datetime.now()

    # instantiate object & variables
    N_sample = int(1e4)
    _A = -10
    _B = 10
    _x_ct = int(1e4)
    n_bins = 50

    h2_arg = [2,1] # [x0,z] for h2 .. the 1/x**2 function

    # do calculation
    print("\n ... get RJM samples with logstic distribution")
    function_g = MyFunctions.pdf_g
    function_h = MyFunctions.logistic_distr
    sample1, ar1 = RejectionMethod.getSamplesRejectionMethod(N_samples=N_sample, func_g=function_g , func_h=function_h , c = 3,
                                              A=_A, B=_B, x_ct=_x_ct, args_h=None,printit=True, plotit=True)
    print("\n ... plotting histograms of accepted values with h1..")
    plt.hist(sample1, bins=n_bins, density=True)
    plt.title("Histogramm of sample1")
    plt.grid()
    plt.show()

    print("\n ... get RJM samples with 1/x**2 envelope")
    function_g = MyFunctions.pdf_g
    function_h = MyFunctions.envelope_func2
    sample2, ar2 = RejectionMethod.getSamplesRejectionMethod(N_samples=N_sample,func_g=function_g ,func_h=function_h, c = 5,
                                              A=_A, B=_B, x_ct=_x_ct,args_h=h2_arg,printit=True, plotit=True)

    print("\n ... plotting histograms of accepted values with h2..")
    plt.hist(sample2, bins=n_bins, density=True)
    plt.title("Histogramm of sample2")
    plt.grid()
    plt.show()

    # END
    print("acceptance rates: ",ar1," ",ar2)
    t2 = datetime.now()
    deltaT = t2-t1
    print("\n-- execution time in seconds: ",int(deltaT.seconds))
    print("\n-- FINISHED --")

