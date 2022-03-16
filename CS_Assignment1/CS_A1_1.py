
# AUTHOR: Christoph Muehleder 01604413
# DATE: March 11, 2021
# Computer Simulations, Assignment 1, Exercise 1

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

class Histograms:

    # Constructor
    def __init__(self):
        # seed random
        np.random.seed(datetime.now().microsecond)

    # A.18
    @staticmethod
    def sigma_hi_bayes(bi, sig_pi):
        return sig_pi/bi

    # A.15
    @staticmethod
    def sigma_pi_bayes(N,nb,pi_exp):
        sig_pi_helper = []
        if pi_exp is int or pi_exp is float:
            return np.sqrt(pi_exp*(1-pi_exp)/(N+nb+2))
        else:
            for pi_expi in pi_exp:
                sig_pi_helper.append(np.sqrt(pi_expi*(1-pi_expi)/(N+nb+2)))
        return np.array(sig_pi_helper)

    # A.14
    @staticmethod
    def pi_exp_bayes(Ni,nb):
        pi_helper = []
        if Ni is int or Ni is int:
            return (Ni+1)/(Ni+1+nb)
        else:
            for ni in Ni:
                pi_helper.append((ni+1)/(ni+1+nb))
        return np.array(pi_helper)

    # A.3
    @staticmethod
    def sigma_Ni(Ni,N):
        t1 = Ni
        t2 = (1-(Ni/N))
        return np.sqrt(t1*t2)

    # A.4
    @staticmethod
    def sigma_Ni_simple(Ni):
        return np.sqrt(Ni)

    # A.5
    @staticmethod
    def sigma_pi_simple(Ni,N):
        return Histograms.sigma_Ni(Ni,N)/N

    @staticmethod
    def getNormedHistogram(N_sample, n_bins=10):

        sample = np.random.random(N_sample)
        hist, bins = np.histogram(sample, n_bins, density=False)
        norm = np.sum(hist)
        hist_normed = hist / norm          # Norm = N * b

        return hist_normed, bins, hist


if __name__ == '__main__':

    print("\n--- executing 1) ---")
    np.random.seed(datetime.now().microsecond)                      # seed random
    t1 = datetime.now()                                             # measure execution time
    histo = Histograms()                                            # instantiate object of this class
    bins_ct = 10                                                       # how many bins
    bins_arr = np.linspace(start=0.0, stop=1.0, num=bins_ct + 1)
    L = 1000                                                        # for frequentist analysis
    Ns = np.array([1e3, 1e4, 1e5],dtype=int)                        # array of samplesizes
    bar = 5                                                         # index of bar to look at
    i = 0                                                           # counter value

    print("--- a) ---")

    pi1, b, Ni = histo.getNormedHistogram(N_sample=Ns[0])            # a) with Histogram plot

    # plot histograms
    plt.hist(np.random.rand(Ns[0]), bins_ct, density=False)
    plt.title("a1) Histogram absolute; N = "+str(Ns[0]))
    plt.ylabel("Ni")
    plt.xlabel("x")
    plt.grid()
    plt.show()
    plt.hist(np.random.rand(Ns[0]), bins_ct, density=False, weights=[1 / (Ns[0]) for i in range(Ns[0])])
    plt.ylabel("pi")
    plt.xlabel("x")
    plt.title("a2) Histogram relative; N = " + str(Ns[0]))
    plt.grid()
    plt.show()

    print("--- b) ---")

    flucs = []                                      # contains L heights of bar#5 probabilites / heights
    hist_widths = []                                # histogram widths (of fluctuations of each L runs of N)

    histograms = []

    print(" ... L = ",L," runs for each samplesize N")
    # make L histograms for each N
    for n in Ns:
        ps1 = []
        for j in range(L):                                      # make a 1000 times
            rands = np.random.rand(n)                           # N randoms between 0 and 1
            pi2, b = np.histogram(rands, bins=bins_ct)          # make histogram of randoms
            pi2 = pi2 / n                                       # norm the histogram
            ps1.append(pi2[bar])                                # append bar height to list

        flucs.append(ps1)

    # calculate histogram widths
    print("... calculating histogram widths")
    for i in range(len(flucs)):
        histo_t, bins_t = np.histogram(a=flucs[i], bins=bins_ct, density=False)
        plt.hist(x=flucs[i], bins =bins_ct, density=False)
        plt.title("Histogram of Fluctuations of Bar #5\nN = "+str(Ns[i]))
        plt.grid()
        plt.show()
        hist_widths.append(bins_t[len(bins_t)-2]-bins_t[1])


    # plot characteristic of bar width against N
    print(" ... plotting histo-width vs N")
    plt.figure()
    plt.title("Fluctuation-Histogram-Width as function of N samples")
    plt.ylabel("W(N)")
    plt.xlabel("N")
    plt.plot(Ns, hist_widths)
    plt.grid()
    plt.show()

    # c) plot histogram with errorbars, calc using frequentist expression
    print("--- c) ---")
    # error bars source https://stackoverflow.com/questions/11774822/matplotlib-histogram-with-errorbars

    sig_pis = []
    pis = []                                # container for every stdev and pi and the bins
    bins = []

    print(" ... calculate error bars")
    # histograms with each N
    for i in range(len(Ns)):
        sample1 = np.random.random(Ns[i])                                        # randoms
        pi, bin = np.histogram(a=sample1, bins=bins_arr)                         # pi of each bin
        pi = pi/Ns[i]

        sig_pis.append(Histograms.sigma_pi_simple(Ni=pi*Ns[i],N=Ns[i]))           # calculate sigma pi for errorbars
        pis.append(pi)
        bins.append(bin)

    # plot everything with errorbars
    for i in range(len(Ns)):

        bincenters = 0.5 * (bins[i][1:] + bins[i][:-1])
        plt.bar(bincenters, pis[i], color='r', width=0.05)                              # draw histogram
        plt.errorbar(bincenters, pis[i], yerr=sig_pis[i])                               # adds errorbars
        plt.title("Histogram of N = " + str(Ns[i]) + " random sample\n with errorbars")
        plt.grid()
        plt.show()

    # --- END ---
    t2 = datetime.now()
    deltaT = t2 - t1
    print("\n-- execution time in seconds: ", int(deltaT.seconds))
    print("\n-- FINISHED --")