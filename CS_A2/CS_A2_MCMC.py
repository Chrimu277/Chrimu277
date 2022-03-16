# AUTHOR: Christoph Muehleder 01604413
# DATE: March 23, 2021
# Computer Simulations, Assignment 2

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

from scipy import integrate as SPI
from libs.CS_A1_1 import Histograms

# SOURCES
# https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

class MCMC:

    # Constructor
    def __init__(self):
        pass

    # the basic gauß distribution
    @staticmethod
    def norm_distr(x,x0,sig):
        t1 = 1/np.sqrt(2*np.pi)
        t2 = ((x-x0)**2)/(2*sig)
        return t1*np.exp(-t2)

    #the PDF we want to random walk in, draw samples from
    @staticmethod
    def pi_x(x,x0,sig):
        return 0.5 * (MCMC.norm_distr(x,x0,sig)+MCMC.norm_distr(x,-x0,sig))

    # draw proposal xb from neighbourhood of xa
    @staticmethod
    def x_beta(x_alpha,dx) -> float:
        return x_alpha + (np.random.rand()-0.5) * dx

    @staticmethod
    def isAccepted(pi_i, pi_j) -> bool:
        if np.random.rand() < (pi_j/pi_i):
            return True
        else:
            return False

    @staticmethod
    def doChain(N, chi, sig, dx, x0 = 0) -> np.ndarray:
        chain = []                                              # the chain to be returned
        accepted = 0                                            # debug purpose

        # metropolis algorithm
        x = (np.random.rand() - 0.5) * dx if x0 == 0 else x0    # starting state

        for j in range(int(N)):                                 # go N steps in chain
            xb = x + (np.random.rand() - 0.5) * dx

            pi_j = MCMC.pi_x(xb, chi, sig)                      # get the proposal probability qij = qji
            pi_i = MCMC.pi_x(x, chi, sig)

            #DECIDE
            if MCMC.isAccepted(pi_i, pi_j):                     # (1.23) with qij = qji = 1
                accepted += 1                                   # accept
                chain.append(xb)                                # switch to new state
                x = xb
            else:                                               # not accepted
                chain.append(x)                                 # stay at old state

        return np.array(chain)

    # gets array of exponential function
    # returns autocorrelation time
    @staticmethod
    def getAutoCorrTime(arr) -> float:
        avg = np.average(arr)                                       # get avg of noise ( bottom line)
        ampl = arr[0]-avg                                           # get difference to peak

        index = 5                                                    # if method fails, 5 works good
        for i in range(len(arr)):                                    # look for index where Amplitude drops by 3/4
            if(arr[i] < ampl/4):
                index = i
                break

        arr_scope = arr[:index]                                      # cut off rest of array
        arr_log = np.log(arr_scope)                                  # take logarithmic of that part
        p, res = np.polyfit(np.arange(len(arr_log)),arr_log,1)       # returns polynomial coefficients and residuals

        # autocorr time is inverse slope !
        return np.abs(1/p)


    # returns index for when a function settles down (minimum variance ansatz)
    @staticmethod
    def get_relaxed_index(arr: np.ndarray, updown: int = 2):
        newarr = arr.copy()
        varminindex = updown
        varmin : float = 1e5
        for i in range(updown,len(newarr)-updown):
            arr_scope = newarr[i-updown:i+updown]
            var = np.var(arr_scope)
            if var < varmin:
                varmin = var
                varminindex = i

        return varminindex


if __name__ == '__main__':

    _SHOWPLOTS = False
    print("\n--- executing MCMC ---")
    np.random.seed(datetime.now().microsecond)                      # seed random
    t1 = datetime.now()   # measure execution time

    # construct State Space x1..x100
    _n_states = 50
    _xn = np.linspace(-15,15,_n_states)

    # parameters for normal distribution
    _sigs = np.array([1,1,1])                  # stdev
    _chis = np.array([0,2,6])                  # X0
    _dx = np.array([2,3,12])                  # 0.55, 0.54 , 0.21 and good curves for 5,10,25
    _N_steps = np.array([int(1e5),int(1e5),int(1e5)])

    # plot the 3 distribution_functions to show
    _refs = []
    _norms = []
    _dx_ref = 1e-3
    _xref = np.arange(-15,15,_dx_ref)
    plt.figure()
    for i in range(len(_sigs)):
        _refs.append(MCMC.pi_x(_xref, _chis[i], _sigs[i]))
        plt.plot(_xref, _refs[i])
        _norms.append(SPI.cumtrapz(y=_refs[i],x=_xref,dx=_dx_ref).max())
    plt.title("Distribution Functions PI(x)")
    plt.legend(["Chi = 0","Chi = 2", "Chi = 6"])
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.savefig("results/a_distributions.png")
    if _SHOWPLOTS : plt.show()
    plt.close()

    # execute 3 runs of markov chain with metropolis algorithm
    print("\n... running Markov Chains")
    chains = []                                # marcov chains
    for i in range(len(_sigs)):
        chain = (MCMC.doChain(_N_steps[i],_chis[i],_sigs[i],_dx[i]))
        chains.append(chain)

    # make Histograms ... does it look like PI(x) ?
    print("\n... making & plotting Histograms")
    Ni_s = []           # histograms
    pi_s = []           # probabilities
    hi_s = []           # normalized bin heigths
    bin_edges_s = []    # bin edges
    for i in range(len(chains)):
        # plot samples histogram
        plt.figure()
        plt.hist(chains[i], bins=_n_states)
        plt.title("Histogram of Metropolis MCMC; Chain #"+str(i)+"\n[chi sig dx N] = "
            "["+str(_chis[i])+" "+str(_sigs[i])+" "+str(_dx[i])+" "+str(_N_steps[i])+"]")
        plt.savefig("results/a_histo_raw" + str(i) + ".png")
        if _SHOWPLOTS : plt.show()
        plt.close()

        # get all important values
        Ni, bin_edges = np.histogram(chains[i], bins=_n_states)             # Ni : histogram
        pi : list = Ni / np.sum(Ni)                                         # pi : probability of 1 bin
        hi : list = pi / (bin_edges[1] - bin_edges[0])                      # bi: normalized height of 1 bin
        nb : int = len(Ni)                                                  # nb: nr of bins
        N : int = np.sum(Ni)                                                # N: sample number
        bi : float = bin_edges[1]-bin_edges[0]                              # bi: binwidth
        bincenters = 0.5 * (bin_edges[1:] + bin_edges[:-1])                 # centers of bins

        hi_s.append(hi)
        pi_s.append(pi)
        bin_edges_s.append(bin_edges)

        # make errors
        pi_exps = [Histograms.pi_exp_bayes(N, Ni[i], nb) for i in range(len(Ni))]                   # <pi>
        sig_pi_bay = [Histograms.sigma_pi_bayes(N, nb, pi_exps[i]) for i in range(len(pi_exps))]    # sigma p_i
        sig_hi_bay = sig_pi_bay/bi                                                                  # sigma h_i

        # plot errors with reference function
        plt.figure()
        plt.errorbar(bincenters,hi,yerr=sig_hi_bay,ecolor="red",elinewidth=1)
        plt.bar(bincenters,hi,width=bi-0.05)
        plt.plot(_xref,_refs[i],'orange')
        plt.title("Plotting Histogram with errorbars of Chain #"+str(i))
        plt.legend(["Reference Distribution Function", "Line+Errorbars", "Bars"])
        plt.ylabel("hi +- sig_hi")
        plt.xlabel("bins")
        plt.savefig("results/a_histo_errorbars" + str(i) + ".png")
        if _SHOWPLOTS : plt.show()
        plt.close()


    # --- b)    make AVGS and STDEV -- are results reasonable ? YES
    print("\n... making AVG and STDEV")
    avgs = []
    std = []
    for i in range(len(chains)):
        avgs.append(np.average(chains[i]))
        std.append(np.std(chains[i]))

    print("\tAVGs: ",avgs)
    print("\tSTD: ",std)

    # --- c)    plot x(t) chain!
    print("\n... plotting x(t) for each chain")
    for i in range(len(chains)):
        plt.figure()
        plt.scatter(np.arange(len(chains[i])),chains[i],s=0.3)
        plt.title("x(t) of Chain #"+str(i))
        plt.savefig("results/c_x_of_t_chain" + str(i) + ".png")
        if _SHOWPLOTS : plt.show()
        plt.close()

    # --- d) autocorrelation function
    # gives how strongly x_n are correlated to each other
    # is related to covariance
    # random sample --> long correlations are BAD
    act = []                            # auto correlation times
    print("\n... getting autocorrelation times with library function")
    for i in range(len(chains)):
        acf = np.correlate(chains[i],chains[i],"full") # length 2*N
        acf = acf / (np.max(acf))
        m = int(len(acf)/2)
        plotdistance = 50

        # plot ac function
        plt.figure()
        plt.plot(np.arange(plotdistance), acf[m:m+plotdistance])
        plt.title("ACF of Chain #" + str(i))
        plt.ylabel("log( acf( t ) )")
        plt.yscale("log")
        plt.xlabel("i")
        plt.savefig("results/d_autocorr_funct"+str(i)+".png")
        if _SHOWPLOTS: plt.show()
        plt.close()

        # get autocorrelation time
        imax = acf.argmax()                 # index of peak
        act.append(MCMC.getAutoCorrTime(acf[imax:]))

    print("\tautocorr. times: ",act)

    # --- e) binning analysis
    # look vor convergence of sigma_Oik (stdev of averages of bins)
    print("\n... binning analysis & plot")
    acf_times_graphical = []
    var_k_converged = []
    for i in range(len(chains)):
        k : float = 1                 # nr values in each bin
        N : int = len(chains[i])    # total nr of values
        new_chain_blocks = []       # list of arrays ( chain split in k parts , for several k's)
        ks = []                     # list of the k's
        var_k_norm = []             # the y values to be plotted

        while k < N/2:                                                      # for every k
            splitchain = np.array_split(chains[i], int(N/k))
            new_chain_blocks.append(splitchain)      # split chain in N/k parts
            ks.append(k)
            N_Bk = len(splitchain)                   # number of blocks

            # from here: formulas from assignment sheet
            O_ik = []
            for block in splitchain:                # calculate the average of each block
                O_ik.append(np.average(block))
            var_k = np.var(O_ik)                      # get variance of averages
            var_k_norm.append(var_k/N_Bk)                 # sigma of O_ik ( averages of blocks )

            # make k bigger
            k *= 1.25

        # plot stdev of Ob - should converge when k is much larger than autocorrelation time
        relaxindex = MCMC.get_relaxed_index(var_k_norm,updown=3)           # looks for index before oscillatiosn start
        acf_times_graphical.append(0.5 * (var_k_norm[relaxindex] / var_k_norm[0]))  # autocorrelation time (1.71)
        var_k_converged.append(var_k_norm[relaxindex])
        # plot the bin analysis
        plt.plot(ks,var_k_norm)
        plt.scatter(ks,var_k_norm)
        plt.axvline(x=ks[relaxindex])
        plt.title("Binning Analysis of Chain #"+str(i+1))
        plt.xlabel("k")
        plt.ylabel("sig_k2 / N_Bk")
        plt.xscale("log")
        plt.savefig("results/e_binning_"+str(i)+".png")
        if _SHOWPLOTS : plt.show()
        plt.close()

    print("\tfound autocorr. times: ",acf_times_graphical)
    print("\tvar_k converges to: ",var_k_converged)

    # --- END ---
    t2 = datetime.now()
    deltaT = t2 - t1
    print("\n-- execution time in seconds: ", int(deltaT.seconds))
    print("-- FINISHED --")