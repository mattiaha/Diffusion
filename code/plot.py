# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:04:00 2020

@author: matti
"""

import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt("CN1data.txt", skiprows = 1)
data2 = np.loadtxt("CN2data.txt", skiprows = 1)
data1ef = np.loadtxt("Ef1data.txt", skiprows = 1)
data2ef = np.loadtxt("Ef2data.txt", skiprows = 1)
data1eb = np.loadtxt("EB1data.txt", skiprows = 1)
data2eb = np.loadtxt("EB2data.txt", skiprows = 1)
x = data1[:,0]
a1 = data1[:,1]
a2 = data2[:,1]
CN1 = data1[:,2]
CN2 = data2[:,2]
EF1 = data1ef[:,2]
EF2 = data2ef[:,2]
EB1 = data1eb[:,2]
EB2 = data2eb[:,2]

n = len(x)
errorEf1 =  np.zeros(n)
errorEf2 =  np.zeros(n)
errorEb1 =  np.zeros(n)
errorEb2 =  np.zeros(n)
errorCN1 =  np.zeros(n)
errorCN2 =  np.zeros(n)
for i in range(1,n):
    errorEf1[i] = np.log10(abs(EF1[i]-a1[i])/a1[i])
    errorEf2[i] = np.log10(abs(EF2[i]-a2[i])/a2[i])
    errorEb1[i] = np.log10(abs(EB1[i]-a1[i])/a1[i])
    errorEb2[i] = np.log10(abs(EB2[i]-a2[i])/a2[i])
    errorCN1[i] = np.log10(abs(CN1[i]-a1[i])/a1[i])
    errorCN2[i] = np.log10(abs(CN2[i]-a2[i])/a2[i])



plt.figure(1)

plt.plot(x,EF1,'b-', label='EulerForward')
plt.title('u(x,t) at $t = 10000 \Delta t$, Forward Euler')
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.show()


plt.figure(2)
plt.plot(x,EF2,'b-', label='EulerForward')
plt.title('u(x,t) at final time, Forward Euler')
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.show()

plt.figure(3)
plt.plot(x[1:-1],errorEf1[1:-1], 'b-', label = 'EulerForward')
plt.plot(x[1:-1],errorEb1[1:-1], 'g-', label = 'EulerBackward')
plt.plot(x[1:-1],errorCN1[1:-1], 'r-', label = 'CN')
plt.legend()
plt.xlabel('x')
plt.ylabel('$\epsilon$ ')
plt.title('$\epsilon$ at t = 10000 $\Delta t$')
plt.show()

plt.figure(4)
plt.plot(x[1:-1],errorEf2[1:-1], 'b-', label = 'EulerForward')
plt.plot(x[1:-1],errorEb2[1:-1], 'g-', label = 'EulerBackward')
plt.plot(x[1:-1],errorCN2[1:-1], 'r-', label = 'CN')
plt.legend()
plt.xlabel('x')
plt.ylabel('$\epsilon$ ')
plt.title('$\epsilon$ at full time')
plt.show()

