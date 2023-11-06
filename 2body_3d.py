"""
Code to solve the two-body problem in the space.
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

def F(Y):
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
    x, y, z, x1, y1, z1 = Y
    # Distance
    r12 = np.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
    # Force
    Fx = - m*M * (x - x1)/r12**3
    Fy = - m*M * (y - y1)/r12**3
    Fz = - m*M * (z - z1)/r12**3
    # Acceleration
    acc = np.array([Fx/m, Fy/m, Fz/m, -Fx/M, -Fy/M, -Fz/M])

    return acc

#=====================================================================
# Numerical integrator, sympletic of 4-th order
#=====================================================================

def Yoshida4(num_steps, tf, f, init, args=()):
    """
    Integrator with Yoshida method

    Parameters
    ----------
    num_steps : int
        number of point of solution
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

    # Some funny coefficents
    l = 2**(1/3)
    w0 = -l/(2-l)
    w1 = 1/(2-l)
    # Other funny coefficents
    c1 = c4 = w1/2
    c2 = c3 = (w0 + w1)/2
    d1 = d3 = w1
    d2 = w0
    # Time steps
    dt = tf/num_steps
    X = np.zeros((num_steps + 1, len(init))) # To store solution
    t = np.zeros(num_steps + 1)              # Time of simutation

    X[0, :] = init                           # Set initial condition

    for i in range(num_steps):               # Loop for solution
        x0 = X[i, ::2]
        v0 = X[i,1::2]

        x1 = x0 + c1*v0*dt
        v1 = v0 + d1*f(x1, *args)*dt

        x2 = x1 + c2*v1*dt
        v2 = v1 + d2*f(x2, *args)*dt

        x3 = x2 + c3*v2*dt
        v3 = v2 + d3*f(x3, *args)*dt

        X[i + 1, ::2] = x3 + c4*v3*dt
        X[i + 1,1::2] = v3
        t[i + 1] = t[i] + dt

    return X, t

#=====================================================================
# Function to compute energy and angular momentum
#=====================================================================

def Ene(x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2):
    """
    Comupute the energy of the system

    Parameters
    ----------
    x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 : 1darray
        coordinates and velocity of body

    Return
    ------
    Ene : 1darray
        Total energy
    """

    e1 = 0.5*m*(vx1**2 + vy1**2 + vz1**2)                      # Kinetic
    e2 = 0.5*M*(vx2**2 + vy2**2 + vz2**2)                      # Kinetic
    e3 = (M*m)/np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)   # Potential

    Ene = e1 + e2 - e3

    return Ene


def Ang(x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2):
    """
    Comupute the angular momentum of the system

    Parameters
    ----------
    x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 : 1darray
        coordinates and velocity of body

    Return
    ------
    ML : 1darray
        angular momentum
    """
    lm_x = m*(y1*vz1 - z1*vy1)
    lm_y = m*(x1*vz1 - z1*vx1)
    lm_z = m*(x1*vy1 - y1*vx1)

    lM_x = M*(y2*vz2 - z2*vy2)
    lM_y = M*(x2*vz2 - z2*vx2)
    lM_z = M*(x2*vy2 - y2*vx2)

    ML = np.array([lm_x + lM_x, lm_y + lM_y, lm_z + lM_z])

    return ML

#=====================================================================
# Computational parameters and initial condition
#=====================================================================

num_steps = 100000
tf = 2

# First mass
x1_0  = -2
vx1_0 = 0
y1_0  = -2
vy1_0 = 10
z1_0  = 0
vz1_0 = -2
# Second mass
x2_0  = 0
vx2_0 = 0
y2_0  = 0
vy2_0 = 0
z2_0  = 0
vz2_0 = 0

init = np.array([x1_0, vx1_0, y1_0, vy1_0, z1_0, vz1_0,
                 x2_0, vx2_0, y2_0, vy2_0, z2_0, vz2_0])

#=====================================================================
# Solution and plot 
#=====================================================================

Sol, t = Yoshida4(num_steps, tf, F, init)
x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 = Sol.T

# Computation of physical quantities
E = Ene(x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2)
L = Ang(x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2)

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
plt.grid()

ax.set_title('2 body problem', fontsize=20)
ax.set_xlabel('X(t)', fontsize=20)
ax.set_ylabel('Y(t)', fontsize=20)
ax.set_zlabel('Z(t)', fontsize=20)

ax.plot(x1, y1, z1, 'k')
ax.plot(x2, y2, z2, 'r')


plt.figure(2)
plt.title('Energy of the system', fontsize=20)
plt.grid()
plt.plot(t, (E-E[0])/E)
plt.xlabel('t', fontsize=20)
plt.ylabel(r'$\frac{E(t)-E(t_0)}{E(t)}$', fontsize=20)

plt.figure(3)
plt.suptitle('Angular momentum: $L(t)-L(t_0)$', fontsize=20)

plt.subplot(311)
plt.grid()
plt.xlabel('t')
plt.ylabel(r'$\frac{Lx(t)-Lx(t_0)}{Lx(t)}$')
plt.plot(t, (L[0, :] - L[0, 0])/L[0, :])

plt.subplot(312)
plt.grid()
plt.xlabel('t')
plt.ylabel(r'$\frac{Ly(t)-Ly(t_0)}{Ly(t)}$')
plt.plot(t, (L[1, :] - L[1, 0])/L[1, :])

plt.subplot(313)
plt.grid()
plt.xlabel('t')
plt.ylabel(r'$\frac{Lz(t)-Lz(t_0)}{Lz(t)}$')
plt.plot(t, (L[2, :] - L[2, 0])/L[2, :])

plt.show()
