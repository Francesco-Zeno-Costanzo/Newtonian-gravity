import numpy as np
import random as rn
import matplotlib.pyplot as plt
from matplotlib import animation

class Corpo:
    '''Classe usata per creare i corpi
    '''

    def __init__(self, x, y, vx, vy, m=1):
        self.m = m
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

    #aggiornamento posizione e velocità con eulero
    def n_vel(self, fx, fy, dt):
        self.vx += fx*dt
        self.vy += fy*dt

    def n_pos(self, dt):
        self.x += self.vx*dt
        self.y += self.vy*dt

class Sistema:
    '''
    Classe per evoluzione del sistema.
    Viene utilizzata la tecnica del softening  per impedire
    divergenze nella foza, sp è il parametro di softening
    '''

    def __init__(self, corpi, G, sp=0):
        self.corpi = corpi
        self.G = G
        self.sp = sp

    def evolvo(self, dt):
        '''
        chimata ad ogni passo temporale, fa evolvere il sistema
        solo di uno step dt, la forza è calcolata secondo la
        legge di gravitazione universale;
        '''

        for corpo_1 in self.corpi:

            fx = 0.0
            fy = 0.0

            for corpo_2 in self.corpi:
                if corpo_1 != corpo_2:

                    dx = corpo_2.x - corpo_1.x
                    dy = corpo_2.y - corpo_1.y

                    d = np.sqrt(dx**2 + dy**2 + self.sp)

                    fx += self.G * corpo_2.m * dx / d**3
                    fy += self.G * corpo_2.m * dy / d**3

            corpo_1.n_vel(fx, fy, dt)

        for corpo in self.corpi:
            corpo.n_pos(dt)

#numero di corpi, pari per come vengono creati
N = 10
C = []
for n in range(N//2):
    #vengono creati due copri alla volta con velocità uguale e opposta
    #per mantenere l'impulso totale del sistema nullo
    v_x = rn.uniform(-0.5, 0.5)
    v_y = rn.uniform(-0.5, 0.5)
    C.append(Corpo(rn.uniform(-0.5, 0.5), rn.uniform(-0.5, 0.5), v_x, v_y))
    C.append(Corpo(rn.uniform(-0.5, 0.5), rn.uniform(-0.5, 0.5), -v_x, -v_y))

G = 1
'''
from scipy import constants as cst
G = cst.value(u'Newtonian constant of gravitation')
and so on for masses and distances for realistic simulation
'''
dim = 2 #il moto avviene su un piano

sist = Sistema(C, G, 0.01)

dt = 1/10000
T = int(1 / dt)
X=np.zeros((dim, T, N))

for t in range(T):
    sist.evolvo(dt)
    for n, corpo in zip(range(N), sist.corpi):
        X[:, t, n] = corpo.x, corpo.y

fig = plt.figure()
plt.title('Problema a n corpi', fontsize=15)
plt.grid()
plt.xlim(np.min(X[0,:,:])-0.5, np.max(X[0,:,:])+0.5)
plt.ylim(np.min(X[1,:,:])-0.5, np.max(X[1,:,:])+0.5)


dot=np.array([])
for corpo in sist.corpi:
    dot=np.append(dot, plt.plot([], [], 'bo'))


def animate(i):
    for k in range(N):
        dot[k].set_data(X[0, i, k], X[1, i, k])

    return dot

anim = animation.FuncAnimation(fig, animate, frames=range(0, T-1, 10), interval=1, blit=True, repeat=True)
#anim.save('N_body.mp4', fps=160,  extra_args=['-vcodec', 'libx264'])
plt.show()
