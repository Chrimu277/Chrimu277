
# AUTHOR: Christoph Muehleder 01604413
# DATE: March 12, 2021
# Computer Simulations, Assignment 1, Exercise 2

# --- global imports ---
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate as scipy_int
from scipy import interpolate as scipy_ip
from datetime import datetime

# --- local imports ---
from CS_Assignment1.Functions import MyFunctionHelperClass as MyFunctions


# write program that generates random numbers according to logistic distribution   OK implemented
# use inverse transformation method
# find corresponding cumulative distribution function CDF = integral PDF
# h(x) is the PDF ?

# SOURCE CODE sampling from arbirtray distr https://gist.github.com/amarvutha/c2a3ea9d42d238551c694480019a6ce1
# SOURCE scipy_stat.t.pdf() and cdf() https://het.as.utexas.edu/HET/Software/Scipy/generated/scipy.stats.t.html

# class for ITM sampling
class InverseTransformation:

    # empty () constructor
    def __init__(self):
        pass

    # THE INVERSE TRANSFORMATION SAMPLE
    # func_pdf: distribution to draw from
    # N_samples: how many samples to draw
    # A: start, B: stop
    # x_ct: steps for x values
    # func_args: extra arguments besides x for pdf
    @staticmethod
    def getItmSampleFromPdf(func_pdf,A,B,x_ct,N_sample,func_args = None):

        np.random.seed(datetime.now().microsecond)              # seed random
        x_arr = np.linspace(A,B,x_ct)                           # make arrays of x values
        dx = x_arr[1]-x_arr[0] # calculate dx                   # calculate dx for integration purposes

        pdf_arr = func_pdf(x_arr) if (func_args is None) else func_pdf(x_arr,func_args)
        norm = np.max(scipy_int.cumtrapz(pdf_arr,x_arr,dx))              # get norming const
        pdf_arr = pdf_arr/norm                                  # norm the pdf

        cdf_x = scipy_int.cumtrapz(y=pdf_arr,x=x_arr,dx=dx)     # alternative integration - already normed & OK
        inverse_cdf = scipy_ip.interp1d(cdf_x, x_arr[:-1])      # inverse of CDF .. choose xi get through H^-1(x) --> x

        x_rand = np.random.rand(10 * N_sample)                  # make some 10 N random values between 0 and 1
        itm_sample = []                                         # samples of size N_sample to be returned

        # draw samples and fill list
        i=0
        while len(itm_sample) < N_sample:
            i = i + 1
            try:
                itm_sample.append(float(inverse_cdf(x_rand[i])))
            except ValueError:
                # could not sample
                continue

        return itm_sample

# ~ end of Class


if __name__ == '__main__':

    print("\n--- executing 2) ---")

    # instantiate class
    itm = InverseTransformation()

    # instantiate some parameters
    t1 = datetime.now()                             # measure execution time
    samplesize= int(1e5)                            # how many samples to draw
    A = 0                                          # lower limit
    B = 10                                          # upper limit
    x_ct = int(1e4)                                 # how finely to split up array
    x_arr = np.linspace(start=A,stop=B,num=x_ct)    # make array of x-vals
    dx = x_arr[2]-x_arr[1]                          # dx for integration purposes
    n_bins = 50                                     # how many bins in histogram

    # make pdf(x)
    print(" ... making & plotting h(x)")
    pdf = MyFunctions.logistic_distr(x_arr)                     # the pdf(x)
    norm = scipy_int.quad(MyFunctions.logistic_distr, A,B)[0]   # get the norm through integration
    pdf = pdf/norm                                      # norm the pdf

    # plot pdf(x)
    plt.plot(x_arr, pdf)
    plt.title("PDF")
    plt.ylabel("h(x)")
    plt.xlabel("x")
    plt.grid()
    plt.show()

    # make CDF(x)
    print(" ... making & plotting H(x)")
    cdf = scipy_int.cumtrapz(pdf, x_arr, dx=dx)  # integrate pdf = cdf

    # plotting CDF (x)
    plt.plot(x_arr[1:], cdf)
    plt.title("CDF")
    plt.ylabel("H(x)")
    plt.xlabel("x")
    plt.grid()
    plt.show()

    # ! get samples !
    print(" ... getting ",samplesize," samples through ITM")
    result = itm.getItmSampleFromPdf(MyFunctions.logistic_distr,A=A,B=B,x_ct=x_ct, N_sample=samplesize)

    # plot histogram of result
    print(" ... plotting result")
    plt.hist(result,n_bins,density=True)
    plt.title("ITM Sample with "+str(n_bins)+" bins and N = "+str(samplesize))
    plt.grid()
    plt.show()

    # END
    t2 = datetime.now()
    deltaT = t2 - t1
    print("\n-- execution time in seconds: ", int(deltaT.seconds))
    print("\n-- FINISHED --")
