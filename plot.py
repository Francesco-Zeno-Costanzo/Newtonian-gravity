"""
In this code we firstly set parameters for simulatioon.
These parameters will be written in a txt file which is
the input for C code.
The second part we load de data of simulation to make
plot ad namimation of the system.
"""
import numpy as np
import random as rn
import matplotlib.pyplot as plt
from matplotlib import animation

# For better plot
import mplhep
plt.style.use(mplhep.style.CMS)

#===========================================================================
# Computational parameters
#===========================================================================

dt   = 1/20000
T    = int(2/dt)
G    = 1
soft = 0.01
"""
#M    = [1000, 1000, 10]
#X    = [0.5, 0, -0.5, 0, -1.5, 0]
#V    = [0,  20, 0,  -20, 0,   40]
#N    = len(M)
M    = [1000, 1]
X    = [0, 0, -2.5, 0]
V    = [0, 0, 0, 15]
N    = len(M)
"""
N = 100
X = []
V = []
M = [1]*N
np.random.seed(69420)
x = np.linspace(0, 2*np.pi, N)

for n in range(N//2):
    '''
    two bodies are created at a time
    with equal and opposite velocity
    to keep the total momentum of the system zero
    '''
    #v_x, v_y = np.random.uniform(-0.5, 0.5), np.random.uniform(-0.5, 0.5)

    #X.append(np.random.uniform(-0.5, 0.5))
    #X.append(np.random.uniform(-0.5, 0.5))
    #V.append(v_x)
    #V.append(v_y)

    #X.append(np.random.uniform(-0.5, 0.5))
    #X.append(np.random.uniform(-0.5, 0.5))
    #V.append(-v_x)
    #V.append(-v_y)

    X.append(np.random.normal(-0.35, 0.05))
    X.append(np.random.normal(0.0, 0.05))
    V.append(0)
    V.append(5)

    X.append(np.random.normal(0.35, 0.05))
    X.append(np.random.normal(0.0, 0.05))
    V.append(0)
    V.append(-5)

# write on a file that will be read by c code
with open("init.txt", "w") as f:
    f.write(f"{G} {soft} {dt} {N} {T} \n")
    for m, x, y, vx, vy in zip(M, X[0::2], X[1::2], V[0::2], V[1::2]):
        f.write(f"{m} {x} {y} {vx} {vy} \n")
    del m, x, vx, vy

#exit()

#===========================================================================
# Load data
#===========================================================================

T = T + 1 # for initial condition

x, y, vx, vy = np.loadtxt("plot.txt", unpack=True)

x  = np.reshape(x,  (T, N))
y  = np.reshape(y,  (T, N))
vx = np.reshape(vx, (T, N))
vy = np.reshape(vy, (T, N))

E, L = np.loadtxt("EL.txt", unpack=True)

#===========================================================================
# Plot and animation
#===========================================================================

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
plt.plot(t, (L - L[0])/L)
plt.xlabel('t')#, fontsize=20)
plt.ylabel(r'$\frac{L(t)-L(t_0)}{L(t)}$')#, fontsize=20)
#plt.savefig("ang_yosh_simpl.pdf")

fig = plt.figure(2)
plt.grid()
plt.xlim(np.min(x)-0.5, np.max(x)+0.5)
plt.ylim(np.min(y)-0.5, np.max(y)+0.5)
colors = plt.cm.cool(np.linspace(0, 1, N))

dot  = np.array([]) # for the planet
line = np.array([]) # to see the trace

for c in colors:
    dot  = np.append(dot,  plt.plot([], [], 'o', c=c))
    line = np.append(line, plt.plot([], [], '-', c=c))   


def animate(i):
    
    for k in range(N):
        #len_trace = 1000

        # Trace of the trajectory       
        #if i > len_trace:
        #    line[k].set_data(x[i-len_trace:i, k], y[i-len_trace:i, k])
        #else:
        #    line[k].set_data(x[:i, k], y[:i, k])
        
        # Point
        dot[k].set_data(x[i, k], y[i, k])
    
    ALL = [*dot]#, *line]
    return ALL

anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, T, 25),
                               interval=1, blit=True, repeat=True)


plt.title('Binary stars and planet')#, fontsize=20)
plt.xlabel('X(t)')#, fontsize=20)
plt.ylabel('Y(t)')#, fontsize=20)

# Ucomment to save the animation
#anim.save('grav1.mp4', fps=120, extra_args=['-vcodec', 'libx264']) 

plt.show()
