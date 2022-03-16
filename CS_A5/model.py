
# AUTHOR: Christoph Mühleder
# TITLE: Computer Simulations Assignment 5
# DATE of last executable: June 15, 21

# file imports

# global imports
import numpy as np
from typing import List
from datetime import datetime as dt

# local imports

# ATOM Class
class ArgonAtom:

    def __init__(self):
        # speed and position
        self.vx: List[float] = []
        self.vy: List[float] = []
        self.x: List[float] = []
        self.y: List[float] = []



class VerletInt:

    def __init__(self):

        np.random.seed(dt.now().microsecond)

        self.TEmp = 300
        self.kb = 1.38e-23
        self.particles: List[ArgonAtom] = []
        self.sigma = 1
        self.epsilon = 1
        self.m_argon: float = 1.0

        self.Ekins : List[float] = []
        self.Eljs : List[float] = []
        self.Epot : List[float] = []

        self.x0_box = 0
        self.x1_box: float = 1000
        self.y0_box: float = 0
        self.y1_box: float = 10
        self.dT: float = 0.1
        self.steps: int = 1000
        self.N_particles: int = 3
        self.g: float = 10


        self.y0max = 1000
        self.v0max = 10

    def initSimulation(self):

        #init particles
        for i in range(self.N_particles):

            # make atom
            atom = ArgonAtom()

            # random x0 & v0
            x0 = np.random.rand()*(self.x1_box-self.x0_box)
            y0 = np.random.rand()*(self.y1_box-self.y0_box)
            vx0 = np.random.rand()*self.v0max
            vy0 = np.random.rand()*self.v0max

            # append first x and y
            atom.x.append(x0)
            atom.y.append(y0)

            # append first vx and vy
            atom.vx.append(vx0)
            atom.vy.append(vy0)

            # add to list
            self.particles.append(atom)

    def start(self,measure_energies = True):

        # first, calc v1/2
        for p in range(len(self.particles)):
            # old velocities
            v0x = self.particles[p].vx[-1]
            v0y = self.particles[p].vy[-1]

            # calculate forces
            total_force = [0.0, 0.0]
            ljf = self.calcLennardJonesForce(p)
            gf = self.calcGravitationalForce()
            total_force[0] += ljf[0] + gf[0]
            total_force[1] += ljf[1] + gf[1]

            # new velocities
            v1x = v0x + 0.5 * self.dT * total_force[0] / self.m_argon
            v1y = v0y + 0.5 * self.dT * total_force[1] / self.m_argon

            # append velocity
            self.particles[p].vx.append(v1x)
            self.particles[p].vy.append(v1y)

        for i in range(self.steps):

            # for every particle
            for p in range(len(self.particles)):

                # position calculation
                part = self.particles[p]

                # v limit
                v_tot = np.sqrt(part.vx[-1] ** 2 + part.vy[-1] ** 2)
                if v_tot > 20 * self.v0max:
                    part.vx[-1] *= 20 / v_tot
                    part.vy[-1] *= 20 / v_tot

                # new position
                x_new = part.x[-1] + self.dT * part.vx[-1]
                y_new = part.y[-1] + self.dT * part.vy[-1]

                # boundaries ...
                if x_new > self.x1_box:
                    part.vx[-1] *= -1
                    x_new =  self.x1_box

                if x_new < self.x0_box:
                    x_new = self.x0_box
                    part.vx[-1] *= -1

                if y_new < self.y0_box:
                    part.vy[-1] *= -1
                    y_new = self.y0_box

                # append new pos
                self.particles[p].x.append(x_new)
                self.particles[p].y.append(y_new)

            # for every particle
            for p in range(len(self.particles)):

                # force calculation
                total_force = [0,0]
                ljf = self.calcLennardJonesForce(p)
                gf = self.calcGravitationalForce()

                total_force[0] += ljf[0] + gf[0]
                total_force[1] += ljf[1] + gf[1]

                # velocity calculation
                part = self.particles[p]
                vx_new = part.vx[-1] + self.dT*total_force[0] / self.m_argon
                vy_new = part.vy[-1] + self.dT*total_force[1] / self.m_argon
                self.particles[p].vx.append(vx_new)
                self.particles[p].vy.append(vy_new)

            if measure_energies:
                self.Eljs.append(self.calcLennardJonesPotAll())
                self.Ekins.append(self.calcKineticEnergyAll())
                self.Epot.append(self.calcGravitationalPotAll())


    def printAtoms(self):

        for i in range(len(self.particles)):
            p = (self.particles[i])
            print(f"{i}: {p.x[-1]}\t\t{p.y[-1]}\t\t"
                  f"{p.vx[-1]}\t\t{p.vy[-1]}")

    # calculate systems kinetic energy
    def calcKineticEnergyAll(self):
        KE = 0.0
        for i in range(len(self.particles)):
            KE += self.calcKineticEnergy(i)
        return KE
    # calculates kinetic energy of particle i
    def calcKineticEnergy(self,i):
        vx = self.particles[i].vx[-1]
        vy = self.particles[i].vy[-1]
        v_tot_sqr = vx*vx+vy*vy

        return 0.5*self.m_argon*v_tot_sqr
        # calculates the system's LJ potential energy

    # calculate system's lennard jones potential energy
    def calcLennardJonesPotAll(self):
        LJ = 0.0
        for i in range(len(self.particles)):
            LJ += self.calcLennardJonesPot(i)
        return LJ
    # calculates the lennard jones potential energy for particle i
    def calcLennardJonesPot(self, i):
        sum = 0.0
        for j in range(len(self.particles)):
            if i == j:
                pass
            else:
                rij= self.calcRij(i,j)[0]
                sum += (self.sigma/rij)**12 - (self.sigma/rij)**6

        return 4*self.epsilon*sum

    # calculates the lennard jones force for particle i
    def calcLennardJonesForce(self, i):
        force_x = 0.0
        force_y = 0.0
        for j in range(len(self.particles)):
            if i == j:
                pass
            else:
                rij,rx,ry = self.calcRij(i, j)

                force_x += (48*rx)/(rij**2) * ((self.sigma / rij) ** 12 - (self.sigma / rij) ** 6)
                force_y += (48*ry)/(rij**2) * ((self.sigma / rij) ** 12 - (self.sigma / rij) ** 6)


        return [-force_x,-force_y]

    # calculate system's gravitational energy
    def calcGravitationalPotAll(self):
        GR = 0.0
        for i in range(len(self.particles)):
            GR += self.calcGravitationalPot(i)
        return GR
    # returns the gravitational potential energy for particle i
    def calcGravitationalPot(self, i):
        return self.m_argon*self.g*self.particles[i].y[-1]

    # returns the gravitational potential energy for particle i
    def calcGravitationalForce(self) -> tuple:
        force_x = 0
        force_y = - self.m_argon * self.g
        return [force_x,force_y]

    # returns distance of particles with index i and j
    def calcRij(self,i,j):
        x1 = self.particles[i].x[-1]
        x2 = self.particles[j].x[-1]

        y1 = self.particles[i].y[-1]
        y2 = self.particles[j].y[-1]

        dx = x2-x1
        dy = y2-y1

        absolute = np.sqrt(dx*dx+dy*dy)
        return [absolute, dx,dy]
