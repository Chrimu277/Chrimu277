
# Author: Christoph Mühleder 01604413
# Date: May 04 2021

# global imports
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

# local imports
from model import TravellingSalesman as TS

def part1():
    # ---- single run ----

    T0 = 1  # Starting Temperature
    q = 1  # Tk = T0 * k^-q

    N = 50  # number of Points
    L = int(5*1e3)  # Measurements per Temperature L ~= N^2


    ts1 = TS()  # make model
    ts1.initRand(N=N, size=1)  # Random Coordinates USE THIS OR OTHER

    print(f" 1) Single run - with {N} random coordinates !")
    best_conf, best_D = ts1.doSimulation(T0=T0, q=q, L=L,
        folder="results/res1",plot_configs_T=True,plot_minD_T=True,plot_errors=True,plot_timeseries=True)  # run sim

    print("...minimum Distance found: ", best_D)
    ts1.plotMap(best_conf, f"Best Configuration Found:\n"
                                 f"D = {np.round(best_D, 5)}",
                "results/res1/bestconfig_found.png")
    # --- END OF  pt2  ----
    # --- END OF pt1 ---

def part2():
    # ---- multiple runs ----

    # start up
    ts2 = TS()
    ts2.initFile("city-positions.txt")

    # SIMULATION parameters
    N_coords: int = len(ts2.c_arr)
    qs = [0.7,1.00,1.5]
    Ls = [N_coords**2,3*N_coords**2,10*N_coords**2]
    T_start = 10

    # save results here
    D_min = np.zeros((len(qs),len(Ls)), dtype=float)
    D_min_all = 10000
    T_run = np.zeros((len(qs), len(Ls)), dtype=float)
    best_configs = []
    best_config_all = np.zeros((len(qs),len(Ls)))


    print(" 2) Multiple runs with RANDOM Positions & different L, q, T_start")

    tstart = dt.datetime.now()
    for i in range(len(qs)):
        best_configs.append([])
        for j in range(len(Ls)):
            ts2.c_arr = best_config_all
            conf, minD = ts2.doSimulation(T0=T_start, q=qs[i], L=Ls[j],
                folder="results/res2",plot_configs_T=False,plot_minD_T=True,plot_errors=True,plot_timeseries=True)

            # add the best configs up
            best_configs[i].append(conf)
            D_min[i][j] = minD

            # update the best config of all simulations
            if minD <= D_min_all:
                D_min_all = minD
                best_config_all = conf

            # measure runtime
            T_run[i][j] = (dt.datetime.now()-tstart).total_seconds()
            tstart = dt.datetime.now()

    # transpose
    D_min2 = D_min.copy().transpose()
    T_run2= T_run.copy().transpose()

    # plot plot plot
    legends_q = [f"q={qs[i]}" for i in range(len(qs))]
    legends_L = [f"L={Ls[i]}" for i in range(len(Ls))]

    # plot Minimum Energies
    print("... plotting minimum Energy of each run")

    for i in range(len(Ls)):
        plt.plot(qs,D_min2[i])
    plt.title("min(E(q)) for Different L")
    plt.xlabel("q")
    plt.ylabel("D_min")
    plt.legend(legends_L)
    plt.savefig("results/res2/min_EofQ_all_simu.png")
    plt.close()
    for i in range(len(qs)):
        plt.plot(Ls, D_min[i])
    plt.title("min(E(L)) for Different q")
    plt.xlabel("L")
    plt.ylabel("min(E(L))")
    plt.legend(legends_q)
    plt.savefig("results/res2/min_EofL_all_simu.png")
    plt.close()

    # plot Runtime of each run
    print("... plotting runtimes")
    for i in range(len(Ls)):
        plt.plot(qs, T_run2[i])
    plt.title("Runtimes")
    plt.xlabel("q")
    plt.ylabel("t_exe")
    plt.legend(legends_L)
    plt.savefig("results/res2/runtimes_q.png")
    plt.close()
    for i in range(len(qs)):
        plt.plot(Ls, T_run[i])
    plt.title("Runtimes")
    plt.xlabel("L")
    plt.ylabel("t_exe")
    plt.legend(legends_q)
    plt.savefig("results/res2/runtimes_L.png")
    plt.close()

    print("...minimum Distance found: ",D_min_all)
    ts2.plotMap(best_config_all,f"Best Configuration Found for\n"
                                f"q={qs}, L={Ls}, T_start={T_start}, D = {np.round(D_min_all,5)}",
                "results/res2/bestconfig_found.png")
    # --- END OF  pt2  ----

if __name__ == '__main__':

    # START
    print("! executing CS A4: Travelling Salesman !\n\n")
    t1 = dt.datetime.now()

    # calc stuff
    #part1()
    part2()

    #makeVid("results/res1","mov1",5)

    # END
    print("\n\n.. reached End of Program")
    print(f".. Execution Time: {(dt.datetime.now()-t1).total_seconds()} seconds")

    # notify when done
    import winsound
    winsound.Beep(666, 2500)  # f[Hz], t[ms]

