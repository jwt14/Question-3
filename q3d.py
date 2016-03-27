#   HPC Q3d
#   Created by Jan Witold Tomaszewski CID: 00833865

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from math import pi,sin,exp,sqrt

L = 1.                                   #Any of these parameters can be changed by user
N_x = range(30,10,-2)
T = 1
N_t = 500
#N_t = range(1000,10000,1000)
d_t = np.linspace(0.00001, 0.001, num=10)
theta = 0
                                     #Be careful! Integers could yield 0
#d_x = L/float(N_x)
alpha = 0.1
    
RMSE_const_t= []
RMSE_const_x = []
    
for i in range(0,10):
    data = subprocess.Popen("./q3 {} {} {} {} {} {}".format(L, N_x[0], N_t*d_t[i], N_t, alpha, theta),
                            stdout=subprocess.PIPE).communicate()[0]
    a=data.split('\n')
    d_x = L/float(N_x[0])
    heat_vector = []
    heat_analytic = []
    error = 0.
    for j in range(0, N_x[0]+1):
        heat_vector.append(eval(a[int( (0.6*N_t-1)*(N_x[0]+2) +1 +j)]))
        heat_analytic.append(sin(pi*j*d_x/L)*exp(-alpha*pi**2*(0.6*d_t[i]*N_t)/L/L))
    print(i)
    for k in range(0, N_x[0]+1):
        error += (heat_analytic[k]-heat_vector[k])**2
    
    RMSE_const_x.append(sqrt(error/(N_x[0]+1)))
    
for i in range(0,9):
    data = subprocess.Popen("./q3 {} {} {} {} {} {}".format(L, N_x[i], T, N_t, alpha, theta),
                            stdout=subprocess.PIPE).communicate()[0]
    a=data.split('\n')
    d_x = L/float(N_x[i])
    
    heat_vector = []
    heat_analytic = []
    error = 0.
    for j in range(0, N_x[i]+1):
        heat_vector.append(eval(a[int( (0.1*N_t-1)*(N_x[i]+2) +1 +j)]))
        heat_analytic.append(sin(pi*j*d_x/L)*exp(-alpha*pi**2*0.1/L/L))
    
    for k in range(0, N_x[i]+1):
        error += (heat_analytic[k]-heat_vector[k])**2
    
    RMSE_const_t.append(sqrt(error/(N_x[i]+1)))
d_x = np.linspace(1/30., 1/10., num=9)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
axes[0].plot(d_t, RMSE_const_x)
axes[0].set_title("Changing time step")
axes[0].set_xlabel('Time step size')
axes[0].set_ylabel('RMS Error')

axes[1].plot(d_x, RMSE_const_t, 'r')
axes[1].set_title("Changing space step")
axes[1].set_xlabel('Number of space steps')

fig.tight_layout()
fig.savefig("RMS.png")
