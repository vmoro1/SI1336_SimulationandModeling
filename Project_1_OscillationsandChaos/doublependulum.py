from matplotlib import animation
from pylab import *

G = 9.8  # gravitational acceleration

def Ekin(osc):
    """Kinetic energy for the double pendulum"""
    return 1 / (2.0 * osc.m * osc.L * osc.L) * (
            osc.p1 * osc.p1 + 2.0 * osc.p2 * osc.p2 - 2.0 * osc.p1 * osc.p2 * cos(osc.q1 - osc.q2)) / (
                   1 + (sin(osc.q1 - osc.q2)) ** 2)

def Epot(osc):
    """Potential energy for the double pendulum"""
    return osc.m * G * osc.L * (3 - 2 * math.cos(osc.q1) - math.cos(osc.q2))


class Oscillator:
    """Class that holds the parameter and state of a double pendulum"""
    
    def __init__(self):
        self.m = 1
        self.L = 1
        self.t = 0
        self.q1 = 0.32#np.pi / 2#0.1*np.pi#0.32
        self.q2 = 0.11#np.pi / 10#0.15*np.pi#0.11
        
        self.p1 = 0
        self.p2 = 0
        self.E = 15
        self.p2 = np.sqrt(self.p2squaredFromH())
        
        
    def p2squaredFromH(self):
        """Calculates p2 from H"""
        return (self.E - Epot(self)) * (1 + (sin(self.q1 - self.q2)) ** 2) * self.m * self.L * self.L
    
        
class Observables:
    """Class for storing observables for an oscillator"""

    def __init__(self):
        self.time = []    # list to store time
        self.q1list = []  # list to store q1
        self.q2list = []  # list to store q2
        self.p1list = []  # list to store p1
        self.p2list = []  # list to store p2
        self.epot = []  # list to store potential energy
        self.ekin = []  # list to store kinetic energy
        self.etot = []  # list to store total energy
        self.poincare_q1 = []  # list to store q1 for Poincare plot
        self.poincare_p1 = []  # list to store p1 for Poincare plot


def dHdp1(q1, q2, p1, p2, m, L):
    """Derivate of H with respect to p1"""  
    return 1 / (2.0 * m * L * L) * (2*p1 - 2.0 * p2 * cos(q1 - q2))/ (1 + (sin(q1 - q2)) ** 2)


def dHdp2(q1, q2, p1, p2, m, L):
    """Derivate of H with respect to p2"""  
    return 1 / (2.0 * m * L * L) * (4.0 * p2  - 2.0 * p1 * cos(q1 - q2)) / (1 + (sin(q1 - q2)) ** 2)

def dHdq1(q1, q2, p1, p2, m, L):
    """"Derivate of H with respect to q1"""
    return 1 / (2.0 * m * L * L) * (
            -2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) + p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(
        q1 - q2) / (1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * 2.0 * sin(q1)


def dHdq2(q1, q2, p1, p2, m, L):
    """Derivate of H with respect to q2"""
    return 1 / (2.0 * m * L * L) * (
            2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) - p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(q1 - q2) / (
                   1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * sin(q2)


class BaseIntegrator:
    """Base class for different types of intergrators"""
    
    def __init__(self):
        self.dt = 0.01  # time step # lägg i init

    def integrate(self,
                  osc,
                  obs,
                  ):
        """ Perform a single integration step """
        self.timestep(osc, obs)

        """ Append observables to their lists """
        obs.time.append(osc.t)
        obs.q1list.append(osc.q1)
        obs.q2list.append(osc.q2)
        obs.p1list.append(osc.p1)
        obs.p2list.append(osc.p2)
        obs.epot.append(Epot(osc))
        obs.ekin.append(Ekin(osc))
        obs.etot.append(Epot(osc) + Ekin(osc))
        
        if len(obs.q2list) > 1:
            if osc.q2 < 0 and (obs.q2list[-2]) > 0 and osc.p2 > 0:
                obs.poincare_p1.append(osc.p1)
                obs.poincare_q1.append(osc.q1)      
            elif osc.q2 > 0 and (obs.q2list[-2]) < 0 and osc.p2 > 0:
                obs.poincare_p1.append(osc.p1)
                obs.poincare_q1.append(osc.q1)
        
            
    def timestep(self, osc, obs):
        """ Implemented by the child classes """
        pass


class RK4Integrator(BaseIntegrator):
    """Integration according to Runge-Kutta"""
    def timestep(self, osc, obs):
        """simulates one time step with Runge-Kutta."""
        dt = self.dt
        osc.t += dt
        
        # More general way to calculate coefficients
        # prev_coeff = [0, 0, 0 , 0]
        # equations = [dHdp1, dHdq1, dHdp2, dHdq2]
        # coefficients = np.zeros((4,4))
        
        # for i in range(len(prev_coeff)):
        #     for eq in range(len(equations)):
        #         if i < (len(prev_coeff) - 1):
        #             coefficients[i, eq] = ((-1)**eq)*dt*equations[eq](osc.q1+(1/2)*prev_coeff[0], osc.q2+(1/2)*prev_coeff[2], osc.p1+(1/2)*prev_coeff[1], osc.p2+(1/2)*prev_coeff[3], osc.m, osc.L)
        #         else:
        #             coefficients[3, eq] = -1*dt*equations[eq](osc.q1+prev_coeff[0], osc.q2+prev_coeff[2], osc.p1+prev_coeff[1], osc.p2+prev_coeff[3], osc.m, osc.L)
        #     prev_coeff = list(coefficients[i,:])
            
        # coeffs_combination = np.array([1/6, 2/6, 2/6, 1/6])
        # osc.q1 += int(coeffs_combination@coefficients[:,0])
        # osc.p1 += int(coeffs_combination@coefficients[:,1])
        # osc.q2 += int(coeffs_combination@coefficients[:,2])
        # osc.p2 += int(coeffs_combination@coefficients[:,3])
        
        
        k1 = dt*dHdp1(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)
        l1 = -dt*dHdq1(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)
        m1 = dt*dHdp2(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)
        n1 = -dt*dHdq2(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)
        
        k2 = dt*dHdp1(osc.q1 + (1/2)*k1, osc.q2 + (1/2)*m1, osc.p1 + (1/2)*l1, osc.p2 + (1/2)*n1, osc.m, osc.L)
        l2 = -dt*dHdq1(osc.q1 + (1/2)*k1, osc.q2 + (1/2)*m1, osc.p1 + (1/2)*l1, osc.p2 + (1/2)*n1, osc.m, osc.L)
        m2 = dt*dHdp2(osc.q1 + (1/2)*k1, osc.q2 + (1/2)*m1, osc.p1 + (1/2)*l1, osc.p2 + (1/2)*n1, osc.m, osc.L)
        n2 = -dt*dHdq2(osc.q1 + (1/2)*k1, osc.q2 + (1/2)*m1, osc.p1 + (1/2)*l1, osc.p2 + (1/2)*n1, osc.m, osc.L)
        
        k3 = dt*dHdp1(osc.q1 + (1/2)*k2, osc.q2 + (1/2)*m2, osc.p1 + (1/2)*l2, osc.p2 + (1/2)*n2, osc.m, osc.L)
        l3 = -dt*dHdq1(osc.q1 + (1/2)*k2, osc.q2 + (1/2)*m2, osc.p1 + (1/2)*l2, osc.p2 + (1/2)*n2, osc.m, osc.L)
        m3 = dt*dHdp2(osc.q1 + (1/2)*k2, osc.q2 + (1/2)*m2, osc.p1 + (1/2)*l1, osc.p2 + (1/2)*n2, osc.m, osc.L)
        n3 = -dt*dHdq2(osc.q1 + (1/2)*k2, osc.q2 + (1/2)*m2, osc.p1 + (1/2)*l2, osc.p2 + (1/2)*n2, osc.m, osc.L)        
                     
        k4 = dt*dHdp1(osc.q1 + k3, osc.q2 + m3, osc.p1 + l3, osc.p2 + n3, osc.m, osc.L)
        l4 = -dt*dHdq1(osc.q1 + k3, osc.q2 + m3, osc.p1 + l3, osc.p2 + n3, osc.m, osc.L)
        m4 = dt*dHdp2(osc.q1 + k3, osc.q2 + m3, osc.p1 + l3, osc.p2 + n3, osc.m, osc.L)
        n4 = -dt*dHdq2(osc.q1 + k3, osc.q2 + m3, osc.p1 + l3, osc.p2 + n3, osc.m, osc.L)  
        
        osc.q1 += (1/6)*(k1 + 2*k2 + 2*k3 +k4)
        osc.p1 += (1/6)*(l1 + 2*l2 + 2*l3 +l4)
        osc.q2 += (1/6)*(m1 + 2*m2 + 2*m3 +m4)
        osc.p2 += (1/6)*(n1 + 2*n2 + 2*n3 +n4)
        
        # if (abs(osc.q2) < 0.1) and osc.p2 > 0:
        #     obs.poincare_q1.append(osc.q1)
        #     obs.poincare_p1.append(osc.p1)
            

class EulerRichardsonIntegrator(BaseIntegrator):
    """Inegration according to the Euler-Richardson method"""
    def timestep(self, osc, obs):
        """simulates one time step with the Euler-Richardson integrator"""
        pass


def sim(integrator,oscillator=Oscillator(), tmax=30.):
    """Simulates the system without animating it. Returns an Observables 
    object containing all of the stored data. Integratior is any of the objects 
    representing integration, currentyl only Runge-Kutta. 
    oscillator is an Oscillator object. tmax is the simulation time """
    
    obs = Observables()
    numframes = int(tmax / integrator.dt) # Number ofupdates in simulation
    
    for i in range(numframes):
        integrator.integrate(oscillator, obs)
    
    return obs


def poincare_plot(integrator):
    """Produces a Poincare plot. integratos is any of the objects 
    representing integration, currentyl only Runge-Kutta."""
    
    different_q = np.array([[0, 0], [1.1, 0], [0.2, 0.3], [1, 0], [0.5, 0.5]]) # Different innitial conditions
    
    plt.figure()
    plt.xlabel('q1')
    plt.ylabel('p1')
    plt.tight_layout()  
    
    for i in range(len(different_q[:,0])):
        oscillator = Oscillator()
        oscillator.q1 = different_q[i, 0]
        oscillator.q2 = different_q[i, 1]
        observables = sim(integrator, oscillator, 30)
        plt.plot(observables.poincare_q1, observables.poincare_p1, 'ro')  
    plt.show()
        

def animate(framenr, osc, obs, integrator, pendulum_lines, stepsperframe):
    """Animation function which integrates a few steps and return a line for the
    pendulum. osc is an Oscillator object, obs in an Obervables object, 
    integrator is the object representing integrationn."""
    for it in range(stepsperframe):
        integrator.integrate(osc, obs)

    x1 = math.sin(osc.q1);
    y1 = -math.cos(osc.q1);
    x2 = x1 + math.sin(osc.q2);
    y2 = y1 - math.cos(osc.q2)
    pendulum_lines.set_data([0, x1, x2], [0, y1, y2])
    return pendulum_lines,


class Simulation:
    """Simulates and animates the system"""

    def run(self,
            integrator,
            tmax=30.,  # final time
            stepsperframe=5, # how many integration steps between visualising frames
            outfile='energy1.pdf'
            ):
        """Runs the simulation and produces relevant plots"""
        numframes = int(tmax / (stepsperframe * integrator.dt))

        plt.clf()
        fig = plt.figure()
        ax = plt.subplot(xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_lines, = ax.plot([], [], lw=5)

        oscillator = Oscillator()
        obs = Observables()
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(fig, animate,  # init_func=init,
                                       fargs=[oscillator, obs, integrator, pendulum_lines, stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        plt.show()
        plt.waitforbuttonpress(30)

        # a)
        # plt.figure()
        # plt.xlabel('time')
        # plt.ylabel('q1')
        # plt.plot(obs.time, obs.q1list)
        # plt.tight_layout()
        
        # plt.figure()
        # plt.xlabel('time')
        # plt.ylabel('q2')
        # plt.plot(obs.time, obs.q2list)
        # plt.tight_layout()  
        
        # plt.figure()
        # plt.xlabel('time')
        # plt.ylabel('p1')
        # plt.plot(obs.time, obs.p1list)
        # plt.tight_layout()  
        
        # plt.figure()
        # plt.xlabel('time')
        # plt.ylabel('p2')
        # plt.plot(obs.time, obs.p2list)
        # plt.tight_layout()   
        
        # plt.figure()
        # plt.xlabel('time')
        # plt.ylabel('energy')
        # plt.plot(obs.time, obs.epot, obs.time, obs.ekin, obs.time, obs.etot)
        # plt.legend(('Epot', 'Ekin', 'Etot'))
        # plt.tight_layout()  # adapt the plot area tot the text with larger fonts 

        # b)
        # plt.figure()
        # plt.xlabel('q1')
        # plt.ylabel('p1')
        # plt.plot(obs.q1list, obs.p1list)
        # plt.tight_layout()  # adapt the plot area tot the text with larger fonts 

        # plt.figure()
        # plt.xlabel('q2')
        # plt.ylabel('p2')
        # plt.plot(obs.q2list, obs.p2list)
        # plt.tight_layout()  # adapt the plot area tot the text with larger fonts 
        
        # c)
        # plt.figure()
        # plt.xlabel('q1')
        # plt.ylabel('p1')
        # plt.plot(obs.poincare_q1, obs.poincare_p1, 'ro')
        # plt.tight_layout()  # adapt the plot area tot the text with larger fonts
        
        plt.show()
        

# simulation = Simulation()
# simulation.run(integrator=RK4Integrator())
        
# sim(RK4Integrator())
# poincare_plot(integrator=RK4Integrator())        


