
# AUTHOR: Christoph Muehleder 01604413
# DATE: March 15, 2021
# Computer Simulations, Assignment 1, Exercise 4

# --- global imports ---
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# --- local imports ---
from CS_Assignment1.CS_A1_1 import Histograms as Hist
from CS_Assignment1.CS_A1_2 import InverseTransformation as Itm
from CS_Assignment1.CS_A1_3 import RejectionMethod as Rjm
from CS_Assignment1.Functions import MyFunctionHelperClass as MyFunctions


if __name__ == '__main__':
    #START
    print("\n--- executing 4) ---")
    t1 = datetime.now()

    # generate random variables with all 3 generators (itm logstic, rjm 1, rjm 2)
    func_h1 = MyFunctions.logistic_distr
    func_h2 = MyFunctions.envelope_func2
    arg_h2 = [1,1]
    c = (3,5)
    x_cnt = int(1e3)
    N0123 = int(1e3)
    A1 = 0
    B1 = 10
    A2 = -10
    B2 = 10

    sample0 = Itm.getItmSampleFromPdf(func_pdf=MyFunctions.logistic_distr,A=A1,B=B1,N_sample= N0123,x_ct= x_cnt)

    sample1, ar1 = Rjm.getSamplesRejectionMethod(N_samples=N0123, func_g=MyFunctions.pdf_g, func_h=func_h1, c=c[0],
                                              A=A2, B=B2, x_ct=x_cnt, args_h=None, printit=False, plotit=False)
    sample2, ar2 = Rjm.getSamplesRejectionMethod(N_samples=N0123, func_g=MyFunctions.pdf_g, func_h=func_h2, c=c[1],
                                            A=A2, B=B2, x_ct=x_cnt, args_h=arg_h2, printit=False, plotit=False)

    # are acceptance rates reasonable ? YES
    print("acceptance rate 1 * c1 = ",ar1," * ", c[0], " = ", ar1*c[0])
    print("acceptance rate 2 * c2 = ", ar2, " * ", c[1], " = ", ar2 * c[1])

    # perform frequentist analysis -- make errorbars for all methods
    Ns = np.array([int(1e3),int(1e4),int(1e5)])
    nbins = 40                                  # analyze 40 bars

    samples0 = []
    samples1 = []                               # 3 samples of 3 methods
    samples2 = []

    Nis0 = []
    Nis1 = []                                   # Nis of 3 methods
    Nis2 = []

    pis0 = []
    pis1 = []                                   # pi of each sample
    pis2 = []

    bins0 = []
    bins1 = []                                  # bins of each sample
    bins2 = []

    sig_nis0 = []
    sig_nis1 = []                               # sigma Ni of each method and sample
    sig_nis2 = []

    sig_pis0_freq = []
    sig_pis1_freq = []                               # sig_pi of each sample and method
    sig_pis2_freq = []
    sig_pis0_bayes = []

    print("... beginning error analysis")
    # for every samplesize: get the samples
    for N in Ns:

        print("... sampling N = ",N)
        # sample N with every method
        samples0.append(Itm.getItmSampleFromPdf(MyFunctions.logistic_distr, A=A1,B=B1, x_ct=x_cnt, N_sample=N))

        samples1.append(Rjm.getSamplesRejectionMethod(N_samples=N, func_g=MyFunctions.pdf_g, func_h=func_h1, c=c[0],
                                                     A=A2,B=B2, x_ct=x_cnt, args_h=None, printit=False, plotit=False)[0])
        samples2.append(Rjm.getSamplesRejectionMethod(N_samples=N, func_g=MyFunctions.pdf_g, func_h=func_h2, c=c[1],
                                                     A=A2,B=B2, x_ct=x_cnt, args_h=arg_h2, printit=False,
                                                     plotit=False)[0])


    # make histograms and get pi and bins
    for i in range(len(Ns)):

        # get histogram of single sample
        Ni0, b0 = np.histogram(samples0[i], nbins)
        pi0 = Ni0 / Ns[i]
        pis0.append(pi0)
        Nis0.append(Ni0)
        bins0.append(b0)

        Ni1, b1 = np.histogram(samples1[i], nbins)
        pi1 = Ni1 / Ns[i]
        pis1.append(pi1)
        Nis1.append(Ni1)
        bins1.append(b1)

        Ni2, b2 = np.histogram(samples2[i], nbins)
        pi2 = Ni2 / Ns[i]
        pis2.append(pi2)
        Nis2.append(Ni2)
        bins2.append(b2)

    # make the errorbars
    for i in range(len(Ns)):
        sig_pis0_freq.append(Hist.sigma_pi_simple(Ni=pis0[i] * Ns[i], N=Ns[i]))
        sig_pis1_freq.append(Hist.sigma_pi_simple(Ni=pis1[i] * Ns[i], N=Ns[i]))     # frequentist anaysis errors
        sig_pis2_freq.append(Hist.sigma_pi_simple(Ni=pis2[i] * Ns[i], N=Ns[i]))

        # bayesian errors
        #sig_pis0_bayes.append(Hist.sigma_pi_bayes(N=Ns[i], nb=nbins, pi_exp=Hist.pi_exp_bayes(Ni=Nis0, nb=nbins)))

        sig_nis0.append(Hist.sigma_Ni(Ni=pis0[i] * Ns[i], N=Ns[i]))
        sig_nis1.append(Hist.sigma_Ni(Ni=pis1[i] * Ns[i], N=Ns[i]))                 # frequentist analysis errors
        sig_nis2.append(Hist.sigma_Ni(Ni=pis2[i] * Ns[i], N=Ns[i]))


    #print everything with errorbars
    for i in range(len(Ns)):
        bincenters0 = 0.5 * (bins0[i][1:] + bins0[i][:-1])
        bincenters1 = 0.5 * (bins1[i][1:] + bins1[i][:-1])
        bincenters2 = 0.5 * (bins2[i][1:] + bins2[i][:-1])

        plt.bar(bincenters0,pis0[i], color='blue', width=(B1-A1)/nbins)
        plt.errorbar(bincenters0, pis0[i], yerr=sig_pis0_freq[i], ecolor='red')
        plt.title("Histogram of Sample 0\n with FREQ. errorbars, N = "+str(Ns[i]))
        plt.ylabel("p(x)")
        plt.xlabel("x")
        plt.show()

        #plt.bar(bincenters0, pis0[i], color='blue', width=(B1 - A1) / nbins)
        #plt.errorbar(bincenters0, pis0[i], yerr=sig_pis0_bayes[i], ecolor='red')
        #plt.title("Histogram of Sample 0\n with BAYES errorbars, N = " + str(Ns[i]))
        #plt.ylabel("p(x)")
        #plt.xlabel("x")
        #plt.show()

        plt.bar(bincenters1, pis1[i], color='blue', width=(B2-A2)/nbins)
        plt.errorbar(bincenters1, pis1[i], yerr=sig_pis1_freq[i], ecolor='red')
        plt.title("Histogram of Sample 1\n with FREQ. errorbars, N = " + str(Ns[i]))
        plt.ylabel("p(x)")
        plt.xlabel("x")
        plt.show()

        plt.bar(bincenters2, pis2[i], color='blue', width=(B2-A2)/nbins)
        plt.errorbar(bincenters2, pis2[i], yerr=sig_pis2_freq[i], ecolor='red')
        plt.title("Histogram of Sample 2\n with FREQ. errorbars, N = " + str(Ns[i]))
        plt.ylabel("p(x)")
        plt.xlabel("x")
        plt.show()

        # plt.bar(bincenters1, Nis0[i], color='blue', width=0.5)
        # plt.errorbar(bincenters1, Nis0[i], yerr=sig_nis0[i], ecolor='red',elinewidth=0.2)
        # plt.title("Histogram of Sample 0\n with errorbars, N = " + str(Ns[i]))
        # plt.ylabel("N(x)")
        # plt.xlabel("x")
        # plt.show()
        #
        # plt.bar(bincenters1, Nis1[i], color='blue', width=0.5)
        # plt.errorbar(bincenters1, Nis1[i], yerr=sig_nis1[i], ecolor='red')
        # plt.title("Histogram of Sample 1\n with errorbars, N = " + str(Ns[i]))
        # plt.ylabel("p(x)")
        # plt.xlabel("x")
        # plt.show()
        #
        # plt.bar(bincenters1, Nis2[i], color='blue', width=0.5)
        # plt.errorbar(bincenters1, Nis2[i], yerr=sig_nis2[i], ecolor='red')
        # plt.title("Histogram of Sample 2\n with errorbars, N = " + str(Ns[i]))
        # plt.ylabel("p(x)")
        # plt.xlabel("x")
        # plt.show()

    t2 = datetime.now()
    deltaT = t2 - t1
    print("\n-- execution time in seconds: ", int(deltaT.seconds))
    print("-- FINISHED --")

