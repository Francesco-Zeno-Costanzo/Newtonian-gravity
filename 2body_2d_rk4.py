"""
Code to solve the two-body problem on the plane.
A symplectic integrator is used and angular momentum
and energy are verified to be conserved.
"""
import numpy as np
import matplotlib.pyplot as plt

m = 1       # first mass  : planet
M = m*1000  # second mass : star

#=====================================================================
# Compute the force of the system
#=====================================================================

def F(t, Y):
    """
    Comute force of system

    Parameters
    ----------
    Y : 1darray
        array of variables

    Return
    ------
    acc : 1darray
        array with the acceleration of the body
    """
    x, vx, y, vy, x1, vx1, y1, vy1 = Y
    # Distance
    r12 = np.sqrt((x - x1)**2 + (y - y1)**2)
    # Force
    Fx = - m*M * (x - x1)/r12**3
    Fy = - m*M * (y - y1)/r12**3
    # equations of motion
    Y_dot = np.array([vx, Fx/m, vy, Fy/m, vx1, -Fx/M, vy1, -Fy/M])

    return Y_dot

#============================================================================
# Itegration: runge-kutta order 4 method
#============================================================================

def RK4(num_steps, t0, tf, f, init, args=()):
    """
    Integrator With Ruge-Kutta 4th order method

    Parameters
    ----------
    num_steps : int
        number of point of solution
    t0 : float
        lower bound of integration
    tf : float
        upper bound of integration
    f : callable
        function to integrate, must accept vectorial input
    init : 1darray
        array of initial condition
    args : tuple, optional
        extra arguments to pass to f

    Return
    ------
    X : array, shape (num_steps + 1, len(init))
        solution of equation
    t : 1darray
        time
    """
    #time steps
    dt = tf/num_steps

    X = np.zeros((num_steps + 1, len(init))) #matrice delle soluzioni
    t = np.zeros(num_steps + 1)              #array dei tempi

    X[0, :] = init                           #condizioni iniziali
    t[0]    = t0

    for i in range(num_steps):
        xk1 = f(t[i], X[i, :], *args)
        xk2 = f(t[i] + dt/2, X[i, :] + xk1*dt/2, *args)
        xk3 = f(t[i] + dt/2, X[i, :] + xk2*dt/2, *args)
        xk4 = f(t[i] + dt, X[i, :] + xk3*dt, *args)
        X[i + 1, :] = X[i, :] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
        t[i + 1] = t[i] + dt

    return X, t

#=====================================================================
# Function to compute energy and angular momentum
#=====================================================================

def Ene(x1, vx1, y1, vy1, x2, vx2, y2, vy2):
    """
    Comupute the energy of the system

    Parameters
    ----------
    x1, vx1, y1, vy1, x2, vx2, y2, vy2 : 1darray
        coordinates and velocity of body

    Return
    ------
    Ene : 1darray
        Total energy
    """

    e1 = 0.5*m*(vx1**2 + vy1**2)                  # Kinetic
    e2 = 0.5*M*(vx2**2 + vy2**2)                  # Kinetic
    e3 = (M*m)/np.sqrt((x1-x2)**2 + (y1-y2)**2)   # Potential

    Ene = e1 + e2 - e3

    return Ene


def Ang(x1, vx1, y1, vy1, x2, vx2, y2, vy2):
    """
    Comupute the angular momentum of the system

    Parameters
    ----------
    x1, vx1, y1, vy1, x2, vx2, y2, vy2 : 1darray
        coordinates and velocity of body

    Return
    ------
    ML : 1darray
        Total angular momentum
    """
    l1 = x1*vy1 - y1*vx1
    l2 = x2*vy2 - y2*vx2
    ML = m*l1 + M*l2

    return ML

#=====================================================================
# Computational parameters and initial condition
#=====================================================================

num_steps = 100000
tf = 0.5

# First mass
x1_0  = -1
vx1_0 = 0
y1_0  = 0
vy1_0 = 20
# Second mass, rest at origin
x2_0  = 0
vx2_0 = 0
y2_0  = 0
vy2_0 = 0

init = np.array([x1_0, vx1_0, y1_0, vy1_0, x2_0, vx2_0, y2_0, vy2_0])

#=====================================================================
# Solution and plot
#=====================================================================

Sol, t = RK4(num_steps, 0, tf, F, init)
x1, vx1, y1, vy1, x2, vx2, y2, vy2 = Sol.T

# Computation of physical quantities
E = Ene(x1, vx1, y1, vy1, x2, vx2, y2, vy2)
L = Ang(x1, vx1, y1, vy1, x2, vx2, y2, vy2)

plt.figure(1)

plt.title('2 body problem', fontsize=20)
plt.xlabel('X(t)', fontsize=20)
plt.ylabel('Y(t)', fontsize=20)

plt.grid()
plt.plot(x1, y1, 'k')
plt.plot(x2, y2, 'r')

plt.figure(2)
plt.title('Energy of the system', fontsize=20)
plt.grid()
plt.plot(t, (E-E[0])/E)
plt.xlabel('t', fontsize=20)
plt.ylabel(r'$\frac{E(t)-E(t_0)}{E(t)}$', fontsize=20)

plt.figure(3)
plt.title('Angular momentum', fontsize=20)
plt.grid()
plt.plot(t, (L -L[0])/L)
plt.xlabel('t', fontsize=20)
plt.ylabel(r'$\frac{L(t)-L(t_0)}{L(t)}$', fontsize=20)

plt.show()
