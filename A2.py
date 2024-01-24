# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:01:34 2021

@author: dkhun
"""

# Import the needed resources
import utils
import numpy as np
import matplotlib.pyplot as plt

# Range of heights in meters and step size
ht_min = 0
ht_max = 30000
ht_del = 100

# Define the heights and calculate the corresponding parameters using the 
# US Standard Atmosphere
height_std = np.arange(ht_min, ht_max, ht_del)
T_std, p_std, rho_std, a_std = utils.atmosphere(height_std)

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 12, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
plt.figure(1)
plt.plot(T, h/1e3)
plt.plot(T_std - 273, height_std/1e3)
plt.show()

for i in range (0, len(p)):
    if p[i] == 925:
        h1 = h[i]
        print(h[i])
    if p[i] == 700:
        h2 = h[i]
        print (h[i])
thickness = h2 - h1
print(thickness)

tv = np.array([])
for j in range(1, len(T)):
    tempK = T[j] + 273
    w_actual = w[j]/1000
    tv = np.append(tv, (tempK * (1 + (0.61 * w_actual))))

utils.mean_temp_calc(Tv, p, 925, 700)
utils.hypsometric_avgT(925, 700, Tv2, Tv2)
    
    