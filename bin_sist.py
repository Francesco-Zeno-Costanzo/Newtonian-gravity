"""
Code to solve the N-body problem on the plane.
We consider a two star system with a planet (tatooine)
A symplectic integrator is used and angular momentum
and energy are verified to be conserved.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# For better plot
import mplhep
plt.style.use(mplhep.style.CMS)

class Body:
    '''
    Calss for generic body
    
    Parameters
    ----------
    x, y : float
        coordinates on plane
    vx, vy : float
        velocity
    m : float, optional, default 1
        mass
    '''
    def __init__(self, x, y, vx, vy, m=1):
        self.m = m     # mass
        self.x = x     # x coordinate
        self.y = y     # y coordinate
        self.vx = vx   # velocity along x
        self.vy = vy   # velocity along y

    # Methods for updating variables
    def n_vel(self, vec):
        self.vx, self.vy = vec

    def n_pos(self, vec):
        self.x, self.y = vec


class System:
    '''
    Class for system evolution.
    The softening technique is used to prevent
    divergences in strength, sp is the softening parameter
    
    Parameters
    ----------
    bodies : list
        list of object from Body class
    G : float
        universal gravitational constant (=1)
    sp : float, optional, default 0
        softening parameter
    '''

    def __init__(self, bodies, G, sp=0):
        self.bodies = bodies # List of alla body
        self.G     = G       # 6.67x10^-11 = 1
        self.sp    = sp      # Softening parameter


    def force(self, X, Y):
        '''
        the force is calculated according to
        the law of universal gravitation
        
        Parameters
        ----------
        X : 1darray
            array of x coordinate of each body
        Y : 1darray
            array of y coordinate of each body
        
        Return        
        ------
        F : 1darray
            force of all bodies on i-th body
        '''
        
        F = []
        N = len(X)
        
        for i in range(N):

            fx = 0.0
            fy = 0.0

            for j, body in enumerate(self.bodies):
                if i != j:

                    dx = X[j] - X[i]
                    dy = Y[j] - Y[i]

                    d = np.sqrt(dx**2 + dy**2 + self.sp)

                    fx += self.G * body.m * dx / d**3
                    fy += self.G * body.m * dy / d**3
            
            F.append([fx, fy])
        
        return np.array(F)
    
    
    def update(self, dt):
        '''
        Function to evolve the system.
        Fourth-order yoshida is used.
        
        Parameter
        ---------
        dt : float
            time step
        '''
        
        # Some funny coefficents
        l  = 2**(1/3)
        w0 = -l/(2-l)
        w1 = 1/(2-l)
        # Other funny coefficents
        c1 = c4 = w1/2
        c2 = c3 = (w0 + w1)/2
        d1 = d3 = w1
        d2 = w0
        
        x0 = np.array([[body.x,  body.y]  for body in self.bodies])
        v0 = np.array([[body.vx, body.vy] for body in self.bodies])
 
        x1 = x0 + c1*v0*dt
        v1 = v0 + d1*self.force(x1[:, 0], x1[:, 1])*dt

        x2 = x1 + c2*v1*dt
        v2 = v1 + d2*self.force(x2[:, 0], x2[:, 1])*dt
        
        x3 = x2 + c3*v2*dt
        v3 = v2 + d3*self.force(x3[:, 0], x3[:, 1])*dt
        
        x4 = x3 + c4*v3*dt
        v4 = v3
        
        for i, body in enumerate(self.bodies):
            body.n_pos(x4[i])
            body.n_vel(v4[i])
 

class Measure:
    '''
    Parameters
    ----------
    bodies : list
        list of object from Body class
    G : float
        universal gravitational constant (=1)
    sp : float, optional, default 0
        softening parameter
    '''
    def __init__(self, bodies, G, sp=0):
        self.bodies = bodies # List of alla body
        self.G     = G       # 6.67x10^-11 = 1
        self.sp    = sp      # Softening parameter
    
    def energy(self):
        ''' Compute the total energy
        '''
        K = 0
        V = 0
        for body in self.bodies:
            K += 0.5*body.m*(body.vx**2 + body.vy**2)
        
        all_body = self.bodies.copy()
        for body_1 in all_body:
            for body_2 in all_body:
                if body_1 != body_2:
                    dx = body_2.x - body_1.x
                    dy = body_2.y - body_1.y
                    d = np.sqrt(dx**2 + dy**2 + self.sp)

                    V += -body_1.m*body_2.m*self.G/d
            all_body.remove(body_1)

        return K + V
    
    def angular(self):
        ''' Compute the total angular momentum
        '''
        L = 0
        for body in self.bodies:
            l_i = body.m*(body.x*body.vy - body.y*body.vx)
            L += l_i
        
        return L

#===========================================================================
# Creating bodies and the system and computational parameters
#===========================================================================

dt = 1/20000
T  = int(2/dt)
E  = np.zeros(T)
L  = np.zeros(T)
G  = 1

# creation of body
C1 = Body(0.5,  0, 0, 20,  int(1e3))
C2 = Body(-0.5, 0, 0, -20, int(1e3))
C3 = Body(-1.5, 0, 0, 40,  int(1e1))
#C4 = Body(1.5, 0, 0, -40,  int(1e1))
C  = [C1, C2, C3]#, C4]
N  = len(C)
X  = np.zeros((2, T, N)) # 2 because the motion is on a plane
V  = np.zeros((2, T, N)) # 2 because the motion is on a plane

# Creation of the system
soft = 0.0
sist = System( C, G, soft)
M    = Measure(C, G, soft)

#===========================================================================
# Evolution
#===========================================================================

for t in range(T):

    L[t] = M.angular() # measure angular momentum
    E[t] = M.energy()  # measure energy

    sist.update(dt)
    for n, body in enumerate(sist.bodies):
        X[:, t, n] = body.x, body.y
        V[:, t, n] = body.vx, body.vy

#===========================================================================
# Plot and animation
#===========================================================================
np.save("datax.npy", X)
np.save("datav.npy", V)

t = np.linspace(0, T*dt, T)
plt.figure(0)#, figsize=(10, 9))
plt.title('Energy of the system')#, fontsize=20)
plt.grid()
plt.plot(t, (E-E[0])/E)
plt.xlabel('t')#, fontsize=20)
plt.ylabel(r'$\frac{E(t)-E(t_0)}{E(t)}$')#, fontsize=20)
#plt.savefig("ene_yosh_simpl.pdf")


plt.figure(1)#, figsize=(10, 9))
plt.title('Angular momentum')#, fontsize=20)
plt.grid()
plt.plot(t, (L -L[0])/L)
plt.xlabel('t')#, fontsize=20)
plt.ylabel(r'$\frac{L(t)-L(t_0)}{L(t)}$')#, fontsize=20)
#plt.savefig("ang_yosh_simpl.pdf")

fig = plt.figure(2)
plt.grid()
plt.xlim(np.min(X[::2, :])-0.5, np.max(X[::2, :])+0.5)
plt.ylim(np.min(X[1::2,:])-0.5, np.max(X[1::2,:])+0.5)
colors = plt.cm.jet(np.linspace(0, 1, N))

dot  = np.array([]) # for the planet
line = np.array([]) # to see the trace

for c in colors:
    dot  = np.append(dot,  plt.plot([], [], 'o', c=c))
    line = np.append(line, plt.plot([], [], '-', c=c))   

def animate(i):
    
    for k in range(N):
        len_trace = 1000

        # Trace of the trajectory       
        if i > len_trace:
            line[k].set_data((X[0, i-len_trace:i, k],), (X[1, i-len_trace:i, k],))
        else:
            line[k].set_data((X[0, :i, k],), (X[1, :i, k],))
        
        # Point
        dot[k].set_data((X[0, i, k],), (X[1, i, k],))
    
    ALL = [*dot, *line]
    return ALL

anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, T, 50),
                               interval=1, blit=True, repeat=True)


plt.title('Binary stars and planet')#, fontsize=20)
plt.xlabel('X(t)')#, fontsize=20)
plt.ylabel('Y(t)')#, fontsize=20)

# Ucomment to save the animation
#anim.save('grav1.mp4', fps=120, extra_args=['-vcodec', 'libx264']) 

plt.show()
