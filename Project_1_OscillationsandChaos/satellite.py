# simple planar pendulum
# Lagrange method:
# H = (1/2) mv^2 + mgh = (1/2) m L^2 (theta')^2 - mgL cos(theta)
# L = (1/2) m L^2 (theta')^2 + mgL cos(theta)
# d/dt dL/d theta' - dL/d theta = 0
# m L^2 theta" = -mgL sin(theta)
# theta" = -(g/L) sin(theta)
#
# Integration by Euler method:
# theta"(t) = -(g/L)*sin(theta(t))
# theta(t+dt) = theta(t)+theta'(t)*dt
# theta'(t+dt) = theta'(t)+theta"(t)*dt
#
# Integration by Euler-Cromer method:
# theta(t+dt) = theta(t)+theta'(t)*dt
# theta"(t+dt) = -(g/L)*sin(theta(t+dt))
# theta'(t+dt) = theta'(t)+theta"(t+dt)*dt
# Advantage: energy drift cancels out over a whole period

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

plt.rcParams['font.size'] = 18

t = 0.	             # start time 
tmax = 10.	         # final time 
dt = 0.002           # time step
x = 1                # initial position 
y = 0                # initial position 
vx = -0.4            # initial velocity
vy = 0.8             # initial velocity
GM = 1               # gravitational constant times planet mass
m = 1                # satelite mass

time = []            # list to store time
posx = []            # list to store x
posy = []            # list to store y
vel = []             # list to store velocity
epot = []            # list to store potential energy
ekin = []            # list to store kinetic energy
etot = []            # list to store total energy

# parameter controling plot interval
stepsperframe = 25
numframes     = int(tmax/(stepsperframe*dt))

fig = plt.figure()
ax  = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(y=0)    # draw a default hline at y=1 that spans the xrange
plt.axvline(x=0)    # draw a default vline at x=1 that spans the yrange
plt.tight_layout()  # adapt the plot area tot the text with larger fonts 
satellite, = ax.plot([], [], 'ro', markersize=20)


def integrate():
    """Perform a single integration step"""
    global t, x, y, vx, vy
    
    # Euler:
    r   = sqrt(x*x + y*y)
    x  += dt*vx
    y  += dt*vy
    fx  = -GM*m*x/(r*r*r) # force along x
    fy  = -GM*m*y/(r*r*r) # force along y
    vx += dt*fx/m
    vy += dt*fy/m

    v2  = vx*vx + vy*vy   # velocity squared
    
    t += dt 
    time.append(t)
    posx.append(x)
    posy.append(y)
    vel.append(sqrt(v2))
    epot.append(  -GM*m/r )
    ekin.append( 0.5*m*v2 )
    etot.append( 0.5*m*v2 - GM*m/r )


def init():
    """Initialize the animation as empty """
    satellite.set_data([], [])
    return satellite,
   
def animate(framenr):
    """Animation function which integrates a few steps and return a line for the pendulum"""
    
    for it in range(stepsperframe):
        integrate()

    xx = (0, x)
    yy = (0, y)
    satellite.set_data(xx, yy)
    return satellite,

# Call the animator, blit=True means only re-draw parts that have changed
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=numframes, interval=50, blit=True, repeat=False)


plt.show()  # show the animation
plt.waitforbuttonpress(timeout=30)

plt.figure()
plt.clf()
plt.xlabel('x')
plt.ylabel('y')
plt.plot(posx,posy)
plt.tight_layout()

plt.figure()
plt.clf()
plt.xlabel('time')
plt.ylabel('velocity')
plt.plot(time,vel)

plt.figure()
plt.xlabel('time')
plt.ylabel('energy')
plt.plot(time,epot,time,ekin,time,etot)
plt.legend( ('Epot','Ekin','Etot') )
plt.tight_layout()
plt.show() # show the three plots
