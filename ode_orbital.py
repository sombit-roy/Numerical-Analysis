from matplotlib import pyplot as plt
import math
import numpy as np

def RK4(N,x0,xN,y10,y20,f1,f2):

    y1 = [None for _ in range(N)]
    y2 = [None for _ in range(N)]
    x = [None for _ in range(N)]
    
    h = (xN-x0) / (N-1)
    y1[0] = y10
    y2[0] = y20
    x[0] = x0

    for n in range(1,N):
            
        k11 = h * f1(x[n-1], y1[n-1], y2[n-1])
        k21 = h * f2(x[n-1], y1[n-1], y2[n-1])
        
        k12 = h * f1(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2)
        k22 = h * f2(x[n-1] + h/2, y1[n-1] + k11/2, y2[n-1] + k21/2)
        
        k13 = h * f1(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2)
        k23 = h * f2(x[n-1] + h/2, y1[n-1] + k12/2, y2[n-1] + k22/2)
    
        k14 = h * f1(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23)
        k24 = h * f2(x[n-1] + h, y1[n-1] + k13, y2[n-1] + k23)
    
        y1[n] = y1[n-1] + (k11 + 2*k12 + 2*k13 + k14)/6
        y2[n] = y2[n-1] + (k21 + 2*k22 + 2*k23 + k24)/6
        x[n] = x[n-1] + h

    return x, y1, y2

l = 2                # azimuthal quantum number
m = 1                # magnetic quantum number
N = 5001             # number of points
epsilon = 0.0001     # error tolerance
a = -0.9999          # x-axis boundary condition 1
b = 0.9999           # x-axis boundary condition 2
ya = (-0.9999)**l    # y-axis boundary condition 1
yb = 0.9999          # y-axis boundary condition 2
alpha = []
alpha.append((yb-ya) / (b-a))
alpha.append(2 * alpha[0])

def f1(x,y,z):
    return z

def f2(x,y,z):
    return (2*x*z - l*(l+1)*y + (m**2)/(1 - x**2)*y) / (1 - x**2)

p = []
x, y, z = RK4(N,a,b,ya,alpha[0],f1,f2)
p.append(y[N-1])
x, y, z = RK4(N,a,b,ya,alpha[1],f1,f2)
p.append(y[N-1])
y[N-1] = p[0]
k = 2

while (abs(y[N-1] - yb) > epsilon):
    alpha.append(alpha[k-2] + (yb-p[k-2]) * (alpha[k-1]-alpha[k-2]) / (p[k-1]-p[k-2]))
    x, y, z = RK4(N,a,b,ya,alpha[k],f1,f2)
    p.append(y[N-1])
    k += 1 

th = np.linspace(0,2*np.pi,N)
costh = np.cos(th)
e = []
for i in costh:
    e.append(np.abs(x-i).argmin())
Plm = []
for j in range(N):
    Plm.append(y[e[j]]**2)

fig = plt.figure()
ax = fig.add_subplot(111,projection='polar')
ax.plot(th,Plm)
ax.set_yticklabels([])
plt.show()