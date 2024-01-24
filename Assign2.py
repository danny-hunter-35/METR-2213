# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 16:02:25 2021

@author: Chase
"""
#Plot 1

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

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
plt.figure(1)
plt.plot(T + 273, h/1e3)
plt.plot(T_std, height_std/1e3)
plt.title('OUN Norman Observations at 00Z 04 July 2020')
plt.xlabel('Temperature (K)')
plt.ylabel('Height MSL (km)')
plt.show()

#Plot 2

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
plt.figure(2)
plt.plot(T + 273, h/1e3)
plt.plot(T_std, height_std/1e3)
plt.title('OUN Norman Observations at 12Z 04 July 2020')
plt.xlabel('Temperature (k)')
plt.ylabel('Height MSL (km)')
plt.show()

#-------------------------------------------------------------------------

#plot 3

# Range of heights in meters and step size
ht_min = 0
ht_max = 30000
ht_del = 100

# Define the heights and calculate the corresponding parameters using the 
# US Standard Atmosphere
height_std = np.arange(ht_min, ht_max, ht_del)
T_std, p_std, rho_std, a_std = utils.atmosphere(height_std)

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
plt.figure(3)
plt.plot(p , h/1e3)
plt.plot(p_std /1e2 , height_std/1e3)
plt.title('OUN Norman Observations at 00Z 04 July 2020')
plt.xlabel('Pressure (hPa)')
plt.ylabel('Height MSL (km)')
plt.show()


#plot 4

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
plt.figure(4)
plt.plot(p , h/1e3)
plt.plot(p_std /1e2 , height_std/1e3)
plt.title('OUN Norman Observations at 12Z 04 July 2020')
plt.xlabel('Pressure (hPa)')
plt.ylabel('Height MSL (km)')
plt.show()

#--------------------------------------------------------------------------

#plot 5

# Range of heights in meters and step size
ht_min = 0
ht_max = 30000
ht_del = 100

# Define the heights and calculate the corresponding parameters using the 
# US Standard Atmosphere
height_std = np.arange(ht_min, ht_max, ht_del)
T_std, p_std, rho_std, a_std = utils.atmosphere(height_std)

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
plt.figure(5)
plt.semilogx(p , h/1e3)
plt.semilogx(p_std /1e2 , height_std/1e3)
plt.title('OUN Norman Observations at 00Z 04 July 2020')
plt.xlabel('Pressure (hPa)')
plt.ylabel('Height MSL (km)')
plt.show()

#plot 6

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
plt.figure(6)
plt.semilogx(p , h/1e3)
plt.semilogx(p_std /1e2 , height_std/1e3)
plt.title('OUN Norman Observations at 12Z 04 July 2020')
plt.xlabel('Pressure (hPa)')
plt.ylabel('Height MSL (km)')
plt.show()
#--------------------------------------------------------------------------



h1 = 0
h2 = 0

for i in range (0, len(p)):
    if p[i] == 925:
        print(h[i], 'm')
        h1 = h[i]
for i in range (0, len(p)):
    if p[i] == 700:
        h2 = h[i]
        print(h[i], 'm')

Thick = h2 - h1
print('This is the thickness of the atmosphere between 925 hPa and 700 hPa at 12z', Thick,'m')

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)

h1 = 0
h2 = 0

for i in range (0, len(p)):
    if p[i] == 925:
        print(h[i], 'm')
        h1 = h[i]
for i in range (0, len(p)):
    if p[i] == 700:
        h2 = h[i]
        print(h[i], 'm')

Thick = h2 - h1
print('This is the thickness of the atmosphere between 925 hPa and 700 hPa at 0z', Thick,'m')

#----------------------------------------------------------------------------------------------

Tv = np.array([]) 

for j in range (1, len(T)):
    tempK = T[j] + 273
    w_actual = w[j]/1000
    Tv = np.append(Tv, (tempK * (1 + (0.61 * w_actual))))

#print(Tv)

Tv2= utils.mean_temp_calc(Tv, p, 925, 700)

print(Tv2)


procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 12, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)

Tv = np.array([]) 

for j in range (1, len(T)):
    tempK = T[j] + 273
    w_actual = w[j]/1000
    Tv = np.append(Tv, (tempK * (1 + (0.61 * w_actual))))

#print(Tv)

Tv2= utils.mean_temp_calc(Tv, p, 925, 700)

print(Tv2)

#----------------------------------------------------------------------------------------------------

procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)



Hypso = utils.hypsometric_avgT(925, 700, 292.81, 292.81)


print(Hypso)


procYear, procMonth, procDay, procHour, procSite = 2020, 7, 4, 12, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)



Hypso2 = utils.hypsometric_avgT(925, 700, 291.89, 291.89)


print(Hypso2)







