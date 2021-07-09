import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

T = 0.1                  # tension of string
mu = 0.1                 # mass density of string
mpt = 0.1                # mass of particle
c = math.sqrt(T/mu)      # speed of wave
L = 30                   # length of string
h = 0.01                 # discrete length step
t = 4751                 # total time
k = 0.01                 # discrete time step
M = int(L/h) + 1
x = np.linspace(-L/2,L/2,M)
a = c*k / h
alpha = (mpt*h) / (T*k*k)

u = [[None for i in range(t)] for j in range(M)]
u = np.array(u)

u[0,:] = 0
u[M-1,:] = 0

p = int((int(L/2) - 1) / h)
q = int((int(L/2) + 1) / h)
pt = int(3*L/(4*h))

u[:,0]=0
u[:,1]=0

for n in range(p+1,q):
    u[n][0] = math.exp(-1 / (1-x[n]**2))                       # initial gaussian wave
    u[n][1] = u[n][0] * (1 + 2*k*c*x[n] / ((1-x[n]**2)**2))    # right moving pulse
    
for n in range(1,t-1): 
    for m in range(1,M-1):
        if (m == pt):
            u[m][n+1] = (u[m+1][n] - 2*(1-alpha)*u[m][n] - alpha*u[m][n-1] + u[m-1][n]) / alpha
        else:
            u[m][n+1] = (a**2)*(u[m+1][n] + u[m-1][n]) + 2*(1-a**2)*u[m][n] - u[m][n-1]            
       
fig = plt.figure()
ax = plt.axes()

def animate(i):
    y = u[:,3*i]
    ax.clear()
    ax.plot(x,y)
    plt.plot([-L/2 + h*pt], [u[pt,3*i]], 'o', c='#1f77b4', markersize=5)
    plt.ylim([-0.4, 0.4])
    plt.xlim([-L/2, L/2])
    
ani = FuncAnimation(fig=fig, func=animate, interval=1, frames=4750//3, repeat=False)
plt.show()