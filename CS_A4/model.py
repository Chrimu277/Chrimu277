# Author: Christoph Mühleder 01604413
# Date: May 04 2021

# global imports
import numpy as np
import datetime as dt
import io
import pandas as pd
import matplotlib.pyplot as plt

# local imports

# class for representing the Salesman's coordinates of cities
class TravellingSalesman:

    def __init__(self):

        # declaring members
        self.c_arr      :np.ndarray =   None
        self.c_list     :list       =   None
        self.D          :np.ndarray =   None
        self.beta       :float =        None
        self.D_avg_k    :np.ndarray =   None
        self.D_avar_k   :np.ndarray =   None

        # seed for random numbers
        np.random.seed(dt.datetime.now().microsecond)

    # initalize model with file
    def initFile(self,path:str):
        self.c_arr = self.readCoordinatesFile(path)                              # initialize coords

    # initialize model with random
    def initRand(self,N:int = 25, size:float = 1):
        self.c_arr = self.makeCoordinatesRandom(N,size)                          # initialize coords

    @staticmethod
    def readCoordinatesFile(path:str) -> list:
        print(f"... trying to read file {path}")
        lines = []
        coordinates = []

        with open(path) as f:
            lines = f.readlines()

        for line in lines:
            pureline = line.rstrip(" ").lstrip(" ").strip("\n")
            tokens = pureline.split(sep="    ")
            c = complex(float(tokens[0]), float(tokens[1]))
            coordinates.append(c)

        print(f"... found {len(coordinates)} coordinates")
        return np.array(coordinates)

    @staticmethod
    def makeCoordinatesRandom(N:int = 25, size:float = 1):
        print("... making random Coordinates")
        coords = []
        for i in range(N):
            coord = complex((np.random.rand()-0.5)*size,(np.random.rand()-0.5)*size)
            coords.append(coord)

        return np.array(coords)

    @staticmethod
    def checkConvergence(arr,eps=0.001,n_last_check=20):
        if(len(arr) > n_last_check):
            if np.abs(np.nanstd(arr[-n_last_check:]) / np.average(arr[-n_last_check:])) <= eps:
                return True
        return False

    # starts simulation
    def doSimulation(self, T0:float, q:float, L:int = 100,
                     folder="results/res1", plot_configs_T=False, plot_errors=False, plot_timeseries=False, plot_minD_T=False
                     ) -> np.ndarray:
        # print
        print(f"--- ! new Simulation with [T0,q,N,L] = [{T0}, {q}, {L}]")

        # most important simulation parameters
        self.T0 = T0
        self.q = q
        self.L = int(L)

        # make arrays - hold error/AVG information for each Temp
        D_avg_k = []
        D_std_k = []

        # init observable measurement, 2D array
        self.D = []

        # -- SIMULATE !!!
        actualDistance = self.calcDistanceFull(self.c_arr.copy())
        opt_config = self.c_arr.copy()      # TEMPORARY total optimum config
        opt_D = 100000                      # TEMPORARY total minimum Distance

        opt_configs = []                    # optimal configs for each Temp
        opt_Ds = []                         # optimal distance for each Temp

        k = 0                               # temp counter
        converged = False
        # make temperatures
        T_k = []

        # for each Temperature
        while not converged:
            T_k.append(self.T0*np.power(float(k),float(-q)))
            print(".. measuring at T = ",np.round(T_k[k],5), f"\t{L} Times..")

            self.D.append(np.zeros((self.L),dtype = float))

            for j in range(int(L)):

                # propose coordinate flip
                dd = self.proposeCoordinateFlip(T=T_k[k],tot_dist_now=actualDistance)

                # get new distance / energy
                actualDistance += dd
                self.D[k][j] = actualDistance    # measure

                # if measured is smaller than global
                if self.D[k][j] <= opt_D:
                    opt_D = self.D[k][j]
                    opt_config = self.c_arr.copy()

            # for this T, get optimums
            opt_configs.append(opt_config)
            opt_Ds.append(opt_D)

            # next Temperature
            k += 1

            # check convergence
            converged = self.checkConvergence(opt_Ds,0.0001,20)

        print("---converged")
        # and plot the Temperature's optimum config
        if plot_configs_T:
            for i in range(k):
            # plot best config at that temperature
                self.plotMap(opt_configs[i],
                             "Travelling Salesman - optimal path found for"
                             f"\nL={self.L}, T={np.round(T_k[i], 5)}, D = {np.round(opt_Ds[i], 5)}",
                             f"{folder}/TRSA_Temp{i}_q{q}_L{L}_opt.png")


        # ANALYSIS

        # calculate errors
        if plot_errors:
            for i in range(k):
                D_avg_k.append(np.average(self.D[i]))
                try:
                    D_std_k.append(np.sqrt(np.average(self.D[i] * self.D[i]) - np.average(self.D[i]) ** 2))
                except:
                    print("ERROR: set std = 1000")
                    D_std_k.append(1000.0)

            # plot averages and std with errorbars - over T
            Tpowminus1 = [1 / T_k[i] for i in range(k)]
            plt.plot(Tpowminus1, D_avg_k)
            plt.errorbar(Tpowminus1,D_avg_k, color='blue', yerr=D_std_k, ecolor='red')
            plt.xlabel("1/T")
            plt.ylabel("<E(t)>(T)")
            plt.title("<E(t)> as function of T; with simple error analysis")
            plt.savefig(folder + f"/AVG+STD_Energy_q{self.q}_L{self.L}.png")
            plt.close()

        # plot full time series - over t
        if plot_timeseries:
            for i in range(k):
                ts = np.arange(self.L) + self.L * i
                plt.plot(ts, self.D[i])
            plt.xlabel("t (1..N*L)")
            plt.ylabel("E(t)")
            plt.title("E(t) Full Time Series\n"
                      f"q = {self.q}; L = {self.L}")
            plt.savefig(folder + f"/TimeseriesFull_q{self.q}_L{self.L}.png")
            plt.close()

        # plotting minimum of E (opt Ds ) over 1/T
        if plot_minD_T:
            Tpowminus1 = [1 / T_k[i] for i in range(k)]
            plt.plot(Tpowminus1, opt_Ds)
            plt.xlabel("1/T")
            plt.ylabel("min(E(T))")
            plt.title(f"Plot min of E(t) for each T\n"
                      f"E_min = {np.min(opt_Ds)}")
            plt.savefig(folder + f"/min_E_q{self.q}_Ts{self.T0}_L{self.L}.png")
            plt.close()

        # return best config and best D
        i_best = np.argmin(opt_Ds)
        return opt_configs[i_best], opt_Ds[i_best]


    # proposes Coordinate Flip between 2 random indices
    # returns: Difference in Energy done
    def proposeCoordinateFlip(self,T:float, tot_dist_now:float = None):

        # propose random indices to flip
        indices = np.random.randint(1, len(self.c_arr),2)  # get random site

        # energy now
        D1 = 0.0

        if tot_dist_now is None:
            D1 = self.calcDistanceFull(self.c_arr)
        else:
            D1 = tot_dist_now

        new_coords = self.reverseArr(arr= self.c_arr, i1= indices[0],i2= indices[1]) # change order
        D2 = self.calcDistanceFull(new_coords) # new energy
        dD = D2-D1 #energy difference

        # get probabilites
        rand = np.random.rand()
        p_flip = np.exp(-dD/T)

        # decision for swap
        if rand < p_flip:
            self.c_arr = new_coords.copy()
        else:
            dD = 0.0

        # return energy difference
        return dD


    # takes: self.c_arr
    # does: opens matplotlib plot to show cities & paths
    @staticmethod
    def plotMap(coords:np.ndarray ,plottitle:str, filename:str):

        Xs = [coords[i].real for i in range(len(coords))]
        Ys = [coords[i].imag for i in range(len(coords))]

        plt.plot(Xs,Ys)
        plt.scatter(Xs, Ys,color='red')
        plt.title(plottitle)

        plt.savefig(filename)
        plt.close()

    @staticmethod
    def calcDistance2Points_complex(c1:complex,c2:complex) -> float:
        print("... calculating Distance of 2 Points (complex input)")
        return np.abs(c1-c2)

    @staticmethod
    def calcDistance2Points_index(arr:np.ndarray,i1:int, i2:int) -> float:
        print("... calculating Distance of 2 Points (index input)")
        try:
            return np.abs(arr[i1]-arr[i2])
        except:
            print("ERROR; returning 0")
            return 0

    @staticmethod
    def reverseArr(arr:np.ndarray,i1:int, i2:int):
        arr_cop = arr.copy()
        arr_cop[i1:i2+1] = np.flip(arr_cop[i1:i2+1])
        return arr_cop

    @staticmethod
    def calcDistanceFull(arr:np.ndarray) -> float:
        arr_cop = arr.copy()

        D = np.abs(arr_cop[-1] - arr_cop[0])
        for i in range(len(arr_cop) - 1):
            D += np.abs(arr_cop[i + 1] - arr_cop[i])
        return D