# AUTHOR: Christoph Mühleder
# TITLE: Computer Simulations Assignment 5
# DATE of last executable: June 15, 21

# imports
import scipy
from model import VerletInt, ArgonAtom
from typing import List
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim


# for animation
def update(i):
    i = int(i)
    __ax.clear()
    for n in range(len(simu.particles)):
        xs = (simu.particles[n].x[:i])
        ys = (simu.particles[n].y[:i])
        __ax.scatter(xs, ys,s=0.5)
    __ax.set_title("Trajectory animation, i="+str(i))
    return __ax,

# maxwellian distribution
def maxwell(v_abs,kTm,A):
    return A*v_abs*np.exp(-v_abs**2/(2*kTm)) / kTm

# barometric height function
def barometric(h,p0,hs):
    return p0*np.exp(-h/hs)

# main

if __name__ == '__main__':

    # for animation
    gif = True
    __Speed = 100



    # --- SIMULATION / EXECUTION ----
    print("making simulation")
    simu = VerletInt()  # make simulation class

    # setting simulation parameters
    print("setting parameters")
    simu.N_particles: int = 5
    simu.x0_box = 0
    simu.x1_box :float= float(simu.N_particles)
    simu.y0_box :float = 0
    simu.y1_box :float = float(simu.N_particles)
    simu.dT : float = 1/1e3            # unit time
    simu.steps : int = int(1e4)

    simu.g : float = 1
    simu.m_argon : float = 1
    simu.v0max : float = 1

    # make particles, plot them, plot 3 energies
    print("initialize particles")
    simu.initSimulation()
    print(simu.calcGravitationalPotAll(),simu.calcLennardJonesPotAll(),simu.calcKineticEnergyAll())

    # run simulation loop -- most time consuming part
    print("running simulation algorithm ...")
    simu.start()

    # ---- ANALYSIS --------
    ktm_results = []

    # PLOT TRAJECTORY
    print("plotting trajectory")
    for i in range(len(simu.particles)):
        plt.scatter(simu.particles[i].x, simu.particles[i].y)
        plt.plot(simu.particles[i].x, simu.particles[i].y)
    plt.grid()
    plt.xlabel("x / sigma")
    plt.ylabel("y / sigma")
    plt.title("Trajectories of particles")
    plt.savefig("results/trajectory.png")
    #plt.show()
    plt.close()

    # ANIMATION / GIF
    animlen = len(simu.particles[0].y)
    if gif:
        print("making animation")
        __fig = plt.figure()
        __ax = plt.axes(xlim=(simu.x0_box,simu.x1_box), ylim=(simu.y0_box,simu.y1_box))
        anima = anim.FuncAnimation(__fig, update, np.arange(1, animlen,__Speed),interval = 1, blit = False)
        anima.save("results/trajectory.gif",fps=100)
        plt.close(__fig)

    # PLOT ENERGIES
    print("plotting total energy")
    TOTAL_ENERGY = [simu.Eljs[i] + simu.Epot[i] + simu.Ekins[i] for i in range(len(simu.Eljs))]
    timeaxis = np.arange(len(simu.Ekins))*simu.dT
    plt.plot(timeaxis, simu.Ekins)
    plt.plot(timeaxis, simu.Epot)
    plt.plot(timeaxis, simu.Eljs)
    plt.plot(timeaxis, TOTAL_ENERGY)
    plt.xlabel("t / unit-time")
    plt.ylabel("E / unit-energy")
    plt.legend(["Ekin","Epot","Elj","E_TOT"])
    plt.savefig("results/all_energies.png")
    #plt.show()
    plt.close()

    # HISTOGRAM OF VELOCITY
    velocities_abs = []
    for p in simu.particles:
        for i in range(len(p.vx)):
            velocities_abs.append(np.sqrt(p.vx[i]**2+p.vy[i]**2))

    h,b = np.histogram(np.array(velocities_abs),2*simu.N_particles)
    midpoints = [(b[i+1]+b[i])/2 for i in range(len(b)-1)]

    # FIT VELOCITY DISTR & PLOT
    from scipy.optimize import curve_fit
    popt, pcov = scipy.optimize.curve_fit(maxwell, midpoints, h, p0=[1,1000])
    print("ktm, A =",popt)
    fitfunc = [maxwell(midpoints[i],*popt) for i in range(len(midpoints))]
    plt.plot(midpoints,h)
    plt.plot(midpoints, fitfunc)
    plt.title(f"Histogram and fit of absolute velocities of all particles\n"
              f"kT/m = {np.round(popt[0],5)}")
    plt.xlabel("|v| / [sigma/unittime]")
    plt.ylabel("#")
    plt.savefig("results/histogram_velocity.png")
    #plt.show()
    plt.close()
    ktm_results.append(popt[0])

    # determine estimate of kbT/m from MEAN VELOCITY
    v_squared = [velocities_abs[i]**2 for i in range(len(velocities_abs))] # = 2kT/m
    avg_v_squared = np.average(v_squared)
    ktm2 = avg_v_squared/2
    ktm_results.append(ktm2)

    # reproduce BAROMETRIC HEIGHT
    heights = []
    for p in simu.particles:
        for i in range(len(p.y)):
            heights.append(p.y[i])

    h, b = np.histogram(np.array(heights), 2 * simu.N_particles)
    midpoints = [(b[i + 1] + b[i]) / 2 for i in range(len(b) - 1)]

    # fit & plot BAROMETRIC HEIGHT
    popt, pcov = scipy.optimize.curve_fit(barometric, midpoints, h,p0=[10000,10])
    ktm3 = popt[1] * simu.g
    ktm_results.append(ktm3)
    print("p0,hs = ", popt, pcov)
    fitfunc = [barometric(midpoints[i], *popt) for i in range(len(midpoints))]
    plt.plot(midpoints, h)
    plt.plot(midpoints, fitfunc)
    plt.title(f"Histogram and fit of barometric eq. of all particles\n"
              f"kT/m = {np.round(ktm3,3)}")
    plt.xlabel("h / sigma")
    plt.ylabel("#")
    plt.savefig("results/histogram_barometric.png")
    #plt.show()


    # save results (estimates of kb T / m )
    np.savetxt("results/kT_m_determined", np.array(ktm_results))

