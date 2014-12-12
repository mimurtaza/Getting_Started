#Phugoid Model
#-------------------

#From Prof Loren Barba's Class

from math import sin, cos, log, ceil
import numpy
import matplotlib.pyplot as plt
#%matplotlib inline  #needed only in the Python Notebooks. Gives the plot in the notebook
                                #instead of opening a new window.
from matplotlib import rcParams
rcParams['font.family']='serif'
rcParams['font.size']=16

#Model parameters

g=9.8           #gravity in m/s2
v_t=30.0      #trim velocity in m/s
C_D=1/40.  #drag coefficient
C_L=1.0      #lift coefficient

#Initial conditions

v0=v_t        #start at the trim velocity
theta0=0.   #initial angle of trajectory
x0=0.0        #horizontal position (arbitrary)
y0=1000.    #initial altitude

# Setting up the equation

def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    Parameters
    ----------------
    u: array of float
       array containing the solution at time n.
    Returns
    -----------
    dudt: array of float
             array containing the RHS given u.
    """
    v=u[0]
    theta=u[1]
    x=u[2]
    y=u[3]
    return numpy.array([-g*sin(theta)-C_D/C_L*g/v_t**2*v**2,-g*cos(theta)/v+g/v_t**2*v,\
                        v*cos(theta),v*sin(theta)])

# Implementing Euler's method

def euler_step(u,f,dt):
    """Returns the solution at the next time-step using Euler's method.
    Parameters
    ----------------
    u: array of float
        solution at the previous time-step.
    f: function
       function to compute the right hand-side of the system of equation.
    dt: float
         time increment.

    Returns:
    -----------
    u_n_plus_1: array of float
                        approximate solution at the next time step.
    """

    return u+dt*f(u)

#

T=100.0                                   #final time
dt=0.1                                      #time increment
N=int(T/dt)+1                          #number of time-steps
t=numpy.linspace(0.,T,N)       #time discretization

# initialize the array containing the solution for each time-step
u = numpy.empty((N, 4))
u[0] = numpy.array([v0, theta0, x0, y0])# fill 1st element with initial values

# time loop - Euler method
for n in range(N-1):
    
    u[n+1] = euler_step(u[n], f, dt)

# get the glider's position with respect to the time
x = u[:,2]
y = u[:,3]

# visualization of the path
plt.figure(figsize=(8,6))
plt.grid(True)
plt.xlabel(r'x', fontsize=18)
plt.ylabel(r'y', fontsize=18)
plt.title('Glider trajectory, flight time = %.2f' % T, fontsize=18)
plt.plot(x,y, 'k-', lw=2);

    
