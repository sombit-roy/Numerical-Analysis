import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

R = 0.04            # radius of ball
K = 0.8             # thermal conductivity of ball
c = 1695            # specific heat of ball
rho = 535           # mass density of ball
M = 101             # no. of discrete radial steps
k = 0.05            # discrete time step
uwater = 100        # temperature of water
uroom = 25          # room temperature
a = K / (c*rho)
h = R / (M-1)
alpha = a*k / (h**2)

x = np.linspace(0,R,M)
u = [[None for i in range(1)] for j in range(M)]
u = np.array(u)
u[:,0] = uroom
u[M-1][0] = uwater

n = 0
while (u[0][n] < 40):
    
    u = np.append(u, np.array([[0 for i in range(1)] for j in range(M)]), axis=1)
    
    for m in range(1,M-1):
        u[m][n+1] = (1-2*alpha)*u[m][n] + alpha*((m-1)*u[m-1][n] + (m+1)*u[m+1][n])/m
    
    u[0][n+1] = u[1][n+1]
    u[M-1][n+1] = uwater
    n += 1

u = u.astype(float)

fig = plt.figure()
ax = plt.axes()

def animate(i):
    rad = np.linspace(0,R,M)
    azm = np.linspace(0, 2*np.pi, M)
    r, th = np.meshgrid(rad, azm)
    z = [u[:,10*i] for _ in range(M)]
    plt.subplot(projection='polar')
    plt.pcolormesh(th, r, z, cmap='inferno', vmin=0, vmax=100, shading='auto')
    plt.yticks([0, 0.25*R, 0.5*R, 0.75*R, R])
    plt.colorbar()
    
ani = FuncAnimation(fig=fig, func=animate, interval=1, frames=3056//10, repeat=False)
plt.show()