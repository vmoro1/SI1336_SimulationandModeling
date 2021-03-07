# Performs a molecular dynamics simulation of particles in 2 dimensions

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import random as rnd


def gaussianRandomNumbers(sigma):
    """Generate two Gaussian random numbers with standard deviation sigma, mean 0"""
    w = 2
    while (w >= 1):
        rx1 = 2 * rnd.random() - 1
        rx2 = 2 * rnd.random() - 1
        w = rx1 * rx1 + rx2 * rx2
    
    w = math.sqrt(-2 * math.log(w) / w);
    return sigma * rx1 * w, sigma * rx2 * w


def thermalize(vx, vy, kineticEnergyPerParticle):
    """Assigns Gaussian distributed velocities given an energy per particle"""
    for i in range(0, len(vx)):
        vx[i], vy[i] = gaussianRandomNumbers(kineticEnergyPerParticle)


def pairEnergy(r):
    """"The pair potential. r is the radius"""
    V = 4 * (((1 / r)**12) - ((1/r)**6))
    return V

             
def pairForce(r):
    """The pair force. r is the radius""" 
    f = 4 * (12*((1/r)**13) - 6*((1/r)**7))
    return f


def pbc_dist(x1, y1, x2, y2, Lx, Ly):
    """Calculate the shortest periodic distance, unit cell [0,Lx],[0,Ly]. 
    Returns the difference along x, along y and the distance. Assumes all 
    particles are within [0,Lx],[0,Ly]"""
    dx = x1 - x2
    dy = y1 - y2
    while dx < -0.5*Lx:
        dx += Lx
    while dx > 0.5*Lx:
        dx -= Lx
    while dy < -0.5*Ly:
        dy += Ly
    while dy > 0.5*Ly:
        dy -= Ly
    return dx, dy, math.sqrt(dx*dx + dy*dy)


n = 24              # Number of particles
mass    = 1.0       # Mass of a particle
invmass = 1.0/mass

numPerRow = 6       # The number of particles per row in the initial grid

# Unit cell size
Lx = numPerRow*1.12
Ly = numPerRow*1.12

kB = 1.0        # Boltzmann constant
T = 0.2         # Temperature
kBT = kB * T   

dt = 0.01       # Time step
nsteps = 20000  # Number if steos

# Use a value of 1 if you want to see in detail how the particles interact
numStepsPerFrame = 100
numFrames = nsteps//numStepsPerFrame



x = []
y = []
vx = []
vy = []
fx = []
fy = []

# Initialize the particle position to a nearly hexagonal lattice
x = []
y = []
for i in range (0,n):
    x.append(Lx*0.95/numPerRow*((i % numPerRow) + 0.5*(i/numPerRow)))
    y.append(Lx*0.95/numPerRow*0.87*(i/numPerRow))

    # Make the velocity lists
    vx.append(0)
    vy.append(0)

    # Make the force lists
    fx.append(0)
    fy.append(0)
    
# Assign velocities from a Gaussian distribution
thermalize(vx, vy, kBT)

fig = plt.figure()
ax  = plt.subplot(xlim=(0, Lx), ylim=(0, Ly))


# Start recording for Ekin and the heat capacity at time 5
startTimeForAveraging = 100
startStepForAveraging = startTimeForAveraging/dt

step = 0

# Initialize observables
sumEpot  = 0
sumEpot2 = 0

outt = []
ekinList = []
epotList = []
etotList = []


def integrate():
    """Perform one MD integration step"""
    global step, n, x, y, v, Ekin, Epot, sumEpot, sumEpot2
    
    # Clear the energy and potential
    Epot = 0
    for i in range(0, n):
        fx[i] = 0
        fy[i] = 0

    # Compute the pair potentials and forces
    for i in range(0,n):
        for j in range(i+1,n):
            dx, dy, r = pbc_dist(x[i],y[i],x[j],y[j],Lx,Ly)
            Epot  += pairEnergy(r)
            fij = pairForce(r)
            fx[i] += fij * dx / r
            fy[i] += fij * dy / r
            fx[j] -= fij * dx / r
            fy[j] -= fij * dy / r

    if step > startStepForAveraging:
        sumEpot  += Epot
        sumEpot2 += Epot * Epot

    Ekin = 0
    for i in range(0,n):

        # At the first step we alread have the "full step" velocity
        if step > 0:
            # Update the velocities with a half step
            vx[i] += fx[i]*invmass*0.5*dt
            vy[i] += fy[i]*invmass*0.5*dt

        Ekin += 0.5*mass*(vx[i]*vx[i] + vy[i]*vy[i])    # Add the kinetic energy of particle i to the total

        # Update the velocities with a half step
        vx[i] += fx[i] *invmass*0.5*dt
        vy[i] += fy[i]*invmass*0.5*dt

        # Update the coordinates
        x[i]  += vx[i] * dt
        y[i]  += vy[i] * dt

        # Put particles back in the unit cell, useful for visualization
        x[i] = x[i] % Lx
        y[i] = y[i] % Ly

    step += 1        # Increase the time step


def integrate_some_steps(framenr):
    """Perform several MD integration steps and record observables"""
    global numStepsPerFrame, start, step, kBT, x, y, v, Ekin, Epot, sv, svv

    for i in range(0, numStepsPerFrame):
        integrate()
    
    t = step*dt
    outt.append(t)
    ekinList.append(Ekin)
    epotList.append(Epot)
    etotList.append(Epot + Ekin)

    if step >= startStepForAveraging and step % 1000 == 0:
        EpotAv  = sumEpot/(step + 1 - startStepForAveraging)
        Epot2Av = sumEpot2/(step + 1 - startStepForAveraging)
        Cv = (Epot2Av - EpotAv * EpotAv) / (kBT * T)
        print('time', t, 'Cv =', Cv)

    return ax.scatter(x, y, s=1500, marker='o', c="r"),


# Call the animator, blit=True means only re-draw parts that have changed
anim = animation.FuncAnimation(fig, integrate_some_steps,
                               frames=numFrames, interval=50, blit=True, repeat=False)

plt.show()
plt.waitforbuttonpress(timeout=20)

plt.figure()
plt.title('Time step = ' + str(dt) + ' and T = ' + str(T))
plt.xlabel('time')
plt.ylabel('energy')
plt.plot(outt, ekinList ,outt, epotList, outt, etotList)
plt.legend( ('Ekin','Epot','Etot') )
plt.show()
