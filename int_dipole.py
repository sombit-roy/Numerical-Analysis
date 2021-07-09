import random
import math
import numpy as np
import matplotlib.pyplot as plt

def dipole(theta,phi,r,R,M):
    
    V = np.zeros(M)
    
    for m in range(M):
        
        X = r[m]*math.sin(theta)*math.cos(phi)
        Y = r[m]*math.sin(theta)*math.sin(phi)
        Z = r[m]*math.cos(theta)
        N=100000
        k=0

        def f(x,y):
            return 1/math.sqrt((x-X)**2 + (y-Y)**2 + Z**2)

        for n in range(N):

            x = -1 + 2*random.random()
            y = -1 + 2*random.random()
            s = math.sqrt(x**2 + y**2)

            if s <= R:
                if y > 0:
                    V[m] = V[m] + f(x,y)
                else:
                    V[m] = V[m] - f(x,y)
                k = k+1
        
        V[m] = math.pi*(R**2)*V[m]/k

    return V

theta = 90
R = 1.0

r1 = 0
r2 = 5
hr = 0.1
r = np.arange(r1, r2+hr, hr)

phi1 = 0
phi2 = 2*math.pi
hphi = math.pi/18
p = np.arange(phi1, phi2+hphi, hphi)

rad, pol = np.meshgrid(r,p)
X = rad*np.cos(pol)
Y = rad*np.sin(pol)
M = int((r2 - r1)/hr + 1)
Z = np.empty((0,int((r2 - r1)/hr + 1)))

for angle in p:
    Z = np.append(Z,[dipole(theta, angle, r, R, M)],axis=0)

fig=plt.figure()
ax=fig.add_subplot(projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
plt.show()