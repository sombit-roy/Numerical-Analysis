import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import math

def RK44(N,x0,xN,y10,y20,y30,y40,f1,f2,f3,f4):
    
    y1 = [[None, None, None] for _ in range(N)]
    y2 = [[None, None, None] for _ in range(N)]
    y3 = [[None, None, None] for _ in range(N)]
    y4 = [[None, None, None] for _ in range(N)]
    x = [None for _ in range(N)]
    
    h = (xN - x0)/(N - 1)
    y1[0] = y10
    y2[0] = y20
    y3[0] = y30 
    y4[0] = y40
    x[0] = x0
    
    for n in range(1,N):
        
        k11 = h * f1(x[n-1], y1[n-1], y2[n-1], y3[n-1], y4[n-1])
        k21 = h * f2(x[n-1], y1[n-1], y2[n-1], y3[n-1], y4[n-1])
        k31 = h * f3(x[n-1], y1[n-1], y2[n-1], y3[n-1], y4[n-1])
        k41 = h * f4(x[n-1], y1[n-1], y2[n-1], y3[n-1], y4[n-1])
        
        k12 = h * f1(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2, y3[n-1] + k31/2, y4[n-1] + h*k41/2)
        k22 = h * f2(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2, y3[n-1] + k31/2, y4[n-1] + h*k41/2)
        k32 = h * f3(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2, y3[n-1] + k31/2, y4[n-1] + h*k41/2)
        k42 = h * f4(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2, y3[n-1] + k31/2, y4[n-1] + h*k41/2)
        
        k13 = h * f1(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2, y3[n-1] + k32/2, y4[n-1] + h*k42/2)
        k23 = h * f2(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2, y3[n-1] + k32/2, y4[n-1] + h*k42/2)
        k33 = h * f3(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2, y3[n-1] + k32/2, y4[n-1] + h*k42/2)
        k43 = h * f4(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2, y3[n-1] + k32/2, y4[n-1] + h*k42/2)
        
        k14 = h * f1(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23, y3[n-1] + k33, y4[n-1] + k43)
        k24 = h * f2(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23, y3[n-1] + k33, y4[n-1] + k43)
        k34 = h * f3(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23, y3[n-1] + k33, y4[n-1] + k43)
        k44 = h * f4(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23, y3[n-1] + k33, y4[n-1] + k43)
        
        y1[n] = y1[n-1] + (k11 + 2*k12 + 2*k13 + k14)/6
        y2[n] = y2[n-1] + (k21 + 2*k22 + 2*k23 + k24)/6
        y3[n] = y3[n-1] + (k31 + 2*k32 + 2*k33 + k34)/6
        y4[n] = y4[n-1] + (k41 + 2*k42 + 2*k43 + k44)/6
        x[n] = x[n-1] + h
    
    return x, y1, y2, y3, y4

def f1(t,y1,y2,y3,y4):
    return y3

def f2(t,y1,y2,y3,y4):
    return y4

def f3(t,y1,y2,y3,y4):
    return -(m2 * (y1-y2)) / (math.sqrt((y1[0]-y2[0])**2 + (y1[1]-y2[1])**2 + (y1[2]-y2[2])**2))**3

def f4(t,y1,y2,y3,y4):
    return -(m1 * (y2-y1)) / (math.sqrt((y1[0]-y2[0])**2 + (y1[1]-y2[1])**2 + (y1[2]-y2[2])**2))**3

m1 = 1.1                                        # mass of first body
m2 = 0.8                                        # mass of second body
r10 = np.array([10,0,0], dtype=float)           # position vector of first body
r20 = np.array([-10,0,0], dtype=float)          # position vector of second body
v10 = np.array([0,0.1,0.1], dtype=float)        # velocity vector of first body
v20 = np.array([0,-0.1,0], dtype=float)         # velocity vector of second body
t0 = 0                                          # start time
tN = 1000                                       # end time
N = 1001                                        # no. of discrete time steps

t, r1, r2, v1, v2 = RK44(N,t0,tN,r10,r20,v10,v20,f1,f2,f3,f4)
r1 = np.array(r1)
r2 = np.array(r2)

fig = plt.figure()
ax = Axes3D(fig=fig, auto_add_to_figure=False)
fig.add_axes(ax)

def animate(i):
    ax.clear()
    ax.set_xlim3d(min(min(r1[:,0]), min(r2[:,0])), max(max(r1[:,0]), max(r2[:,0])))
    ax.set_ylim3d(min(min(r1[:,1]), min(r2[:,1])), max(max(r1[:,1]), max(r2[:,1])))
    ax.set_zlim3d(min(min(r1[:,2]), min(r2[:,2])), max(max(r1[:,2]), max(r2[:,2])))
    ax.plot3D(r1[2*i,0], r1[2*i,1], r1[2*i,2], 'o', markersize=15)
    ax.plot3D(r2[2*i,0], r2[2*i,1], r2[2*i,2], 'o', markersize=5)
    
ani = FuncAnimation(fig=fig, func=animate, interval=1, frames=1000//2, repeat=False)
plt.show()