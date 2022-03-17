import numpy as np
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

class Energia:
    '''
    energia del sistema al tempo t
    '''
    def __init__(self, corpi, G, sp=0):
        self.corpi = corpi
        self.G = G
        self.sp = sp

    def cinetica(self):
        K = 0
        for corpo in self.corpi:
            K += 0.5*corpo.m*(corpo.vx**2 + corpo.vy**2)
        return K

    def potenziale(self):
        V = 0
        corpi=self.corpi.copy()
        for corpo_1 in corpi:
            for corpo_2 in corpi:
                if corpo_1 != corpo_2:
                    dx = corpo_2.x - corpo_1.x
                    dy = corpo_2.y - corpo_1.y
                    d = np.sqrt(dx**2 + dy**2 + self.sp)

                    V += -corpo_1.m*corpo_2.m*self.G/d
            corpi.remove(corpo_1)

        return V


#creazione dei corpi nella simulazione
C1 = Corpo(0.5, 0, 0, 20, int(1e3))
C2 = Corpo(-0.5, 0, 0, -20, int(1e3))
C3 = Corpo(-1.5, 0, 0, 40, int(1e1))

C = [C1, C2, C3]

G = 1
'''
from scipy import constants as cst
G = cst.value(u'Newtonian constant of gravitation')
and so on for masses and distances for realistic simulation
'''
N = len(C)
dim = 2 #il moto avviene su un piano

soft = 0.0
sist = Sistema(C, G, soft)
Ene = Energia(C, G, soft)

dt = 1/10000
T = int(1/dt)
X = np.zeros((dim, T, N))
ene = np.zeros(T)

for t in range(T):

    ene[t] = Ene.cinetica() + Ene.potenziale()

    sist.evolvo(dt)
    for n, corpo in zip(range(N), sist.corpi):
        X[:, t, n] = corpo.x, corpo.y

plt.figure(1)
plt.title('Energia del sistema: $E-E(t_0)$', fontsize=20)
plt.plot(np.linspace(0, T, len(ene)), ene-ene[0])
plt.grid()

fig = plt.figure(2)
plt.title('Sistema binario con pianeta', fontsize=15)
plt.grid()
plt.xlim(np.min(X[0,:,:])-0.5, np.max(X[0,:,:])+0.5)
plt.ylim(np.min(X[1,:,:])-0.5, np.max(X[1,:,:])+0.5)

colors = plt.cm.jet(np.linspace(0, 1, N))

dot=np.array([])
line=np.array([])
for c, corpo in zip(colors, sist.corpi):
    dot=np.append(dot, plt.plot([], [], 'o', c=c))
    line=np.append(line, plt.plot([], [], '-', c=c))


def animate(i):
    for k in range(N):
        l_scia = 500 #linghezza della scia lascita sul grafico
        if i>l_scia:
            line[k].set_data(X[0, i-l_scia:i, k], X[1, i-l_scia:i, k])
        else:
            line[k].set_data(X[0, :i, k], X[1, :i, k])

        dot[k].set_data(X[0, i, k], X[1, i, k])

    return dot[0], dot[1], dot[2], line[0], line[1], line[2]

anim = animation.FuncAnimation(fig, animate, frames=range(0, T-1, 5), interval=10, blit=True, repeat=True)
#anim.save('bin_sist.mp4', fps=160,  extra_args=['-vcodec', 'libx264'])

plt.show()
