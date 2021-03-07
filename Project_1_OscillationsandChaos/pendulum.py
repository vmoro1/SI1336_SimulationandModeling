from matplotlib import animation
from pylab import * # Ändra så att det ser snyggt ut

# Global constants
G = 9.8  # gravitational acceleration


class Oscillator:
    """Class for an general, simple oscillator"""
    
    def __init__(self):
        self.m = 1      # mass of the pendulum bob
        self.c = 4      # c=g/L=4
        self. L = G / self.c  # string length
    
        self.t = 0              # the innitial time
        self.theta = np.pi / 2  # the innitial position/angle. When theta=0 the pendlum hangs from the origin in the negative y-direction
        self.dtheta = 0         # the  innitial velocity


class Observables:
    """Class for storing observables for an oscillator"""
    def __init__(self):
        self.time = []    # list to store time
        self.pos = []     # list to store positions
        self.vel = []     # list to store velocities
        self.energy = []  # list to store energy


class BaseSystem:
    """Type of system. Either a harmonic oscillator or a pendlum"""
    
    def force(self, osc):
        """ Implemented by the childclasses  """
        pass


class Harmonic(BaseSystem):
    """A harmonic oscillator"""
    def force(self, osc):
        """Force on the pendlum for a harmonic oscillator. osc is an oscillator object"""
        return -osc.m*osc.c * osc.theta
        
        # Dämpande kraft
        # gamma = 1
        # omega_0 = 2
        # return -osc.m*((omega_0**2)*osc.theta + gamma*osc.dtheta)


class Pendulum(BaseSystem):
    """A class representing a pendlum. osc is an oscillator object"""
    
    def force(self, osc):
        """Force on the pendlum for a pendlum"""
        return -osc.m*osc.c * np.sin(osc.theta)
    
        # Dämpande kraft
        # gamma = 1
        # return -osc.m*(osc.c * np.sin(osc.theta) + gamma*osc.dtheta)


class BaseIntegrator:
    """Base class for different types of intergrators"""
    def __init__(self):
        self.dt = 0.01  # time step

    def integrate(self, simsystem, osc, obs):
        """ Perform a single integration step """
        self.timestep(simsystem, osc, obs)

        # Append observables to their lists
        obs.time.append(osc.t)
        obs.pos.append(osc.theta)
        obs.vel.append(osc.dtheta)
        obs.energy.append(
            0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * osc.theta ** 2)  # harmonic oscillator        

    def timestep(self, simsystem, osc, obs):
        """ Implemented by the childclasses """
        pass

class EulerCromerIntegrator(BaseIntegrator):
    """Integrator according to EulerCromer"""
    
    def timestep(self, simsystem, osc, obs):
        """Solution for one time step according to EulerCramer"""
        osc.t += self.dt       
        accel = simsystem.force(osc) / osc.m 
        
        osc.dtheta = osc.dtheta - osc.c*osc.theta*self.dt # update theta
        osc.theta = osc.theta + osc.dtheta*self.dt        # update dtheta


class VerletIntegrator(BaseIntegrator):
    """Integrator according to Verlet method"""
    
    def timestep(self, simsystem, osc, obs):
        """Solution for one time step according to the Verlet method"""
        osc.t += self.dt
        accel = simsystem.force(osc) / osc.m
        
        osc.theta = osc.theta + osc.dtheta*self.dt + (1/2)*accel*(self.dt**2)
        new_accel = simsystem.force(osc) / osc.m   # Acceleration for next timestep
        osc.dtheta = osc.dtheta + (1/2)*(new_accel + accel)*self.dt
        

class RK4Integrator(BaseIntegrator):
    """Integration according to Runge-Kutta"""
    def timestep(self, simsystem, osc, obs):
        """Solution for one time step according to the Runge-Kutta method"""
        dt = self.dt
        osc.t += dt

        accel = simsystem.force(osc) / osc.m
        force = self.force_pendulum_RK
        
        # Runge-Kutta coefficients
        a1 = accel*dt
        b1 = osc.dtheta*dt
        a2 = (force(osc.m, osc.c, osc.theta + b1/2) / osc.m)*dt
        b2 = (osc.dtheta + a1/2)*dt
        a3 = (force(osc.m, osc.c, osc.theta + b2/2) / osc.m)*dt
        b3 = (osc.dtheta + a2/2)*dt
        a4 = (force(osc.m, osc.c, osc.theta + b3) / osc.m)*dt
        b4 = (osc.dtheta + a3)*dt
        
        osc.dtheta = osc.dtheta + (1/6)*(a1 + 2*a2 + 2*a3 + a4)
        osc.theta = osc.theta+ (1/6)*(b1 + 2*b2 + 2*b3 + b4)
        
    def force_HO_RK(self, mass, c, theta):
        """Force for harmonic oscillator"""
        return -mass*c*theta
    
    def force_pendulum_RK(self, mass, c, theta):
        """Force pendulum"""
        return -mass * c * np.sin(theta)



def animate(framenr, simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe):
    """Animation function which integrates a few steps and return a line for the pendulum"""
    
    for it in range(stepsperframe):
        integrator.integrate(simsystem, oscillator, obs)

    x = np.array([0, np.sin(oscillator.theta)])
    y = np.array([0, -np.cos(oscillator.theta)])
    pendulum_line.set_data(x, y)
    return pendulum_line,


# 1.1
    
def analytical_HO(theta_0, title):
    """Analytical solution to the hramonic oscillator"""
    
    A = theta_0
    g = 9.8
    c = 2     # c=sqrt(g/L)=2
    l = g / (c**2)
    m = 1
    
    
    time = [0]
    theta = [theta_0]
    theta_dot = [0]
    E = [(1/2)*m*g*l*(theta_0**2)]
    
    dt = 0.01
    tmax = 30
    num_steps = int(tmax/dt)
    t = 0
    
    # Store and update values
    for i in range(num_steps):
        t += dt
        theta.append(A*np.cos(c*t))
        theta_dot.append(-A*c*np.sin(c*t))
        E.append((1/2)*m*(l**2)*(theta_dot[-1]**2) + (1/2)*m*g*l*(theta[-1]**2))
        time.append(t)
        
    # Plot the solution
    plt.figure()
    plt.title(title)
    plt.plot(time, theta, label="Position")
    plt.plot(time, theta_dot, label="Velocity")
    plt.plot(time, E, label="Energy")
    plt.xlabel('time')
    plt.legend()
    plt.show()
    


def sim(simsystem,integrator, oscillator=Oscillator(), tmax=30.):
    """Simulates the system without animating it. Returns an Observables 
    object containing all of the stored data. Simsystem is either a Pendulum or 
    Harmonic object. Integratior is any of the objects representing integration 
    methods. oscillator is an Oscillator object. tmax is the simulation time """
    
    obs = Observables()
    numframes = int(tmax / integrator.dt)
    
    for i in range(numframes):
        integrator.integrate(simsystem, oscillator, obs)
    
    return obs
        
#1.2   
def period(simsystem, integrator, tmax=30, stepsperframe=5):
    """Plots the period for versus the innitial theta. Simsystem is either a Pendulum or 
    Harmonic object. Integratior is any of the objects representing integration 
    methods. oscillator is an Oscillator object. tmax is the simulation time"""
    
    
    theta_0 = 0
    theta_final = np.pi - 0.01
    num_thetas = 10
    delta_theta = (theta_final - theta_0) / num_thetas
    theta = 0
    periods = []
    
    for i in range(num_thetas): # Iterate over all initial conditions
        oscillator = Oscillator()
        oscillator.theta = theta
        obs = sim(simsystem, integrator, oscillator)
        theta += delta_theta
        
        thetas = obs.pos
        time = obs.time
        count = 0
        times = []
        prev_time = 0
        
        for i in range(len(thetas) - 1):
            if thetas[i-1] <= thetas[0] <= thetas[i+1]:
                count += 1          
                times.append(time[i] - prev_time)
                prev_time = times[-1]
        
        
        period = sum(times) / count # Average period
        periods.append(period)
    
    plt.figure()
    plt.xlabel('Theta')
    plt.ylabel('Time')
    plt.plot(thetas, periods)
    plt.show()
    

def damped_pendulum():
    """Phase portrait for the damped pendulum. Theta dot is plotted versus theta"""
    
    gamma = 1
    obs = sim(Pendulum(), VerletIntegrator())
    theta = obs.pos
    dtheta = obs.vel
    plt.plot(theta, dtheta)
    

class Simulation:
    """Simulate and animate the system"""

    def run(self,
            simsystem,
            integrator,
            tmax=30.,  # final time
            stepsperframe=5,  # how many integration steps between visualising frames
            title="simulation",  # Name of output file and title shown at the top
            ): 
        """Runs the simulation and produces relevant plots"""
        
        oscillator = Oscillator()
        obs = Observables()
        numframes = int(tmax / (stepsperframe * integrator.dt)) # Number of steps in simulation

        plt.clf()
        # fig = plt.figure()
        ax = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_line, = ax.plot([], [], lw=5)
        plt.title(title)
        # Call the animator
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        plt.show()
        plt.waitforbuttonpress(30)
        
        plt.clf()
        plt.title(title)
        plt.plot(obs.time, obs.pos, label="Position")
        plt.plot(obs.time, obs.vel, label="Velocity")
        plt.plot(obs.time, obs.energy, label="Energy")
        plt.xlabel('time')
        plt.legend()
        plt.savefig(title + ".pdf")
        plt.show()


# 1.1 ======================================== sätt kraften till icke dämpande
simulation = Simulation()

# simulation.run(simsystem=Harmonic(), integrator=EulerCromerIntegrator(), title="Harmonic-EulerCromer")
# simulation.run(simsystem=Pendulum(), integrator=EulerCromerIntegrator(), title="Pendulum-EulerCromer")

# simulation.run(simsystem=Harmonic(), integrator=VerletIntegrator(), title="Harmonic-Verlet")
# simulation.run(simsystem=Pendulum(), integrator=VerletIntegrator(), title="Pendulum-Verlet")
        
# simulation.run(simsystem=Harmonic(), integrator=RK4Integrator(), title="Harmonic-Runge Kutta")  # Se till så att det är rätt kraft, pendel vs HO
# simulation.run(simsystem=Pendulum(), integrator=RK4Integrator(), title="Pendulum-Runge Kutta")

# analytical_HO(0.5*np.pi, 'Analytical Harmonic')
# analytical_HO(0.1*np.pi, 'Analytical Harmonic')

#=============================================================================
  
# 1.2 ====
    
# per = period(simsystem=Harmonic(), integrator=VerletIntegrator())
# per = period(simsystem=Pendulum(), integrator=VerletIntegrator())
        
# ============================================================================

# 1.3 sätt kraften till dämpande i Harmonicklassen ===========================
        
# observables = sim(simsystem=Harmonic(), integrator=VerletIntegrator())
# plt.title('Damped Harmonic Oscillator, gamma = ')
# plt.plot(observables.time, observables.pos, label="Position")
# plt.plot(observables.time, observables.vel, label="Velocity")
# plt.plot(observables.time, observables.energy, label="Energy")
# plt.xlabel('time')
# plt.legend()
# plt.show()      
 
# simulation = Simulation()       
# simulation.run(simsystem=Harmonic(), integrator=VerletIntegrator(), title="Damped Harmonic Oscillator Verlet")
        
#=============================================================================
        
# 1.4  sätt kraften till dämpande i Pendulumlassen ===========================

# observables = sim(simsystem=Pendulum(), integrator=VerletIntegrator(), tmax=60)
# plt.figure()
# plt.title('Phase space portrait for damped pendulum')
# plt.xlabel('theta')
# plt.ylabel('theta dot')
# plt.plot(observables.pos, observables.vel)
# plt.show()

#=============================================================================



# obs = sim(simsystem=Harmonic(), integrator=VerletIntegrator())

# per = period(simsystem=Harmonic(), integrator=VerletIntegrator())
# per = period(simsystem=Pendulum(), integrator=VerletIntegrator())

