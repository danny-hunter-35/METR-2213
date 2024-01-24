# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:47:34 2021

@author: dkhun
"""

# Import the needed resources
import utils
import numpy as np
import matplotlib.pyplot as plt

# Range of heights in meters and step size
ht_min = 0
ht_max = 20000
ht_del = 200

# Define the heights and calculate the corresponding parameters using the 
# US Standard Atmosphere
height = np.arange(ht_min, ht_max, ht_del)
T, p, rho, a = utils.atmosphere(height)
TC = T - 273.15

# Plot the results and generate output files
plt.plot(TC, height/1000)
plt.ylim(ht_min/1000, ht_max/1000)
plt.xlabel('Temperature (C)')
plt.ylabel('Geopotential Height (km)')
plt.title('U.S. Standard Atmosphere')
plt.savefig('std_atm_fig1.png')
plt.show()

plt.plot(p/100, height/1000)
plt.ylim(ht_min/1000, ht_max/1000)
plt.xlabel('Pressure (hPa)')
plt.ylabel('Geopotential Height (km)')
plt.title('U.S. Standard Atmosphere')
plt.savefig('std_atm_fig2.png')
plt.show()

plt.plot(rho, height/1000)
plt.ylim(ht_min/1000, ht_max/1000)
plt.xlabel('Density (kg m$^{-3}$)')
plt.ylabel('Geopotential Height (km)')
plt.title('U.S. Standard Atmosphere')
plt.savefig('std_atm_fig3.png')
plt.show()

ind1 = np.argmin(np.abs(p - 70000))
p1 = p[ind1]
T1 = T[ind1]
ind2 = np.argmin(np.abs(p - 50000))
p2 = p[ind2]
T2 = T[ind2]
Delta_z_hypso = utils.hypsometric_avgT(p1, p2, T1, T2)
z1 = height[ind1]
z2 = height[ind2]
Delta_z_actual = z2 - z1
print(Delta_z_hypso)
print(Delta_z_actual)
