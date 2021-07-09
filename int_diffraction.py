import math
import numpy as np
import matplotlib.pyplot as plt

def Simpson(F,x1,xN,N):
    
    h = (xN-x1) / (N-1)
    I = (F(x1) + F(xN)) * h/3
    x = [None for _ in range(N)]
    x[0] = x1
    
    for j in range(1,N):
        x[j] = x[j-1] + h

    for i in range(1,N-1):
        if (i%2 == 0):
            I = I + 4*h*F(x[i])/3
        else:
            I = I + 2*h*F(x[i])/3

    return I

def f1(t):
    return math.cos((math.pi/2) * t**2)

def f2(t):
    return math.sin((math.pi/2) * t**2)

t = np.linspace(-15,15,1500)
c = np.array([None for _ in range(1500)])
s = np.array([None for _ in range(1500)])
intensity = np.array([None for _ in range(1500)])

for i in range(1500):
    c[i] = Simpson(f1,0,t[i],10000)
    s[i] = Simpson(f2,0,t[i],10000)
    intensity[i] = (0.5-c[i])**2 + (0.5-s[i])**2

diff = np.tile(intensity, (1500,1))
diff = diff.astype(float)
plt.imshow(diff, cmap='gray')
plt.show()