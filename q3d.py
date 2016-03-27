#   HPC Q2c
#   Created by Jan Witold Tomaszewski CID: 00833865

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from math import pi,sin,exp,sqrt

L = 1.                                                                          #Any of these parameters can be changed by user
N_x = range(30,10,-2)
T = 1
N_t = 500
d_t = np.linspace(0.00001, 0.001, num=10)
alpha = 0.1
theta = 0.5                                                                     #Added theta parameter
RMSE_const_t= []
RMSE_const_x = []
    
for i in range(0,len(d_t)):                                                     #First loop: varying time steps
    data = subprocess.Popen("./q3 {} {} {} {} {} {}".format(L, N_x[0], N_t*d_t[i], N_t, alpha, theta),
                            stdout=subprocess.PIPE).communicate()[0]            #Loading data using subprocess routine
    a=data.split('\n')
    d_x = L/float(N_x[0])
    heat_vector = []
    heat_analytic = []
    error = 0.
    for j in range(0, N_x[0]+1):
        heat_vector.append(eval(a[int( (0.6*N_t-1)*(N_x[0]+2) +1 +j)]))
        heat_analytic.append(sin(pi*j*d_x/L)*exp(-alpha*pi**2*(0.6*d_t[i]*N_t)/L/L))

    for k in range(0, N_x[0]+1):
        error += (heat_analytic[k]-heat_vector[k])**2                           #Accummulating the errors
    
    RMSE_const_x.append(sqrt(error/(N_x[0]+1)))                                 #Storing RMS error
    
for i in range(0,len(N_x)):                                                     #Second loop: varying space steps
    data = subprocess.Popen("./q3 {} {} {} {} {} {}".format(L, N_x[i], T, N_t, alpha, theta),
                            stdout=subprocess.PIPE).communicate()[0]            #Loading data using subprocess routine
    a=data.split('\n')
    d_x = L/float(N_x[i])
    
    heat_vector = []
    heat_analytic = []
    error2 = 0.
    for j in range(0, N_x[i]+1):
        heat_vector.append(eval(a[int( (0.1*N_t-1)*(N_x[i]+2) +1 +j)]))
        heat_analytic.append(sin(pi*j*d_x/L)*exp(-alpha*pi**2*0.1/L/L))
    
    for k in range(0, N_x[i]+1):
        error2 += (heat_analytic[k]-heat_vector[k])**2                          #Accummulating the errors
    
    RMSE_const_t.append(sqrt(error2/(N_x[i]+1)))                                #Storing RMS error
d_x = np.linspace(1./N_x[0], 1./N_x[-1], num=len(d_t))


#Plotting routines
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
axes[0].plot(d_t, RMSE_const_x)
axes[0].set_title("Changing time step")
axes[0].set_xlabel('Time step size')
axes[0].set_ylabel('RMS Error')
axes[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
axes[1].plot(d_x, RMSE_const_t, 'r')
axes[1].set_title("Changing space step")
axes[1].set_xlabel('Space step size')
axes[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig.tight_layout()
fig.savefig("RMS_CN.pdf")