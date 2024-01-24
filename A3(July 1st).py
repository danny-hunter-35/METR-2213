import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units
import matplotlib.pyplot as plt
import numpy as np
import utils

procYear, procMonth, procDay, procHour, procSite = 2011, 7, 1, 0, 72357
soundingInfo, p, h, T, Td, RH, w, wdir, sknt, thta, thte, thtv = \
    utils.getWyoming(procYear, procMonth, procDay, procHour, procSite)
p = p[1:]
h = h[1:]
T = T[1:]
Td = Td[1:]
sknt = sknt[1:]
wdir = wdir[1:]

# Setting up the skew t plot
fig = plt.figure(figsize=(9, 9))
skew = SkewT(fig, rotation=30)
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
ww = np.array([0.0002, 0.0004, 0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.045, 0.060]).reshape(-1, 1)
pp = np.linspace(100, 1000) * units.mbar
skew.plot_mixing_lines(ww, pp)
skew.plot(p*units.hPa, T*units.degC, 'r', linewidth=1.5)
skew.plot(p*units.hPa, Td*units.degC, 'g', linewidth=1.5)
plt.xlabel("Temperature (°C)")
plt.ylabel("Pressure (hPa)")
plt.title(soundingInfo)

# Parcel path (black line) and LCL point
parcel_prof = mpcalc.parcel_profile(p*units.hPa, T[0]*units.degC, \
Td[0]*units.degC).to('degC')
lcl_pressure, lcl_temperature = mpcalc.lcl(p[0]*units.hPa, \
T[0]*units.degC, Td[0]*units.degC)
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
skew.plot(p*units.hPa, parcel_prof, 'k', linewidth=2)

# Plot 10 black dots (parcel path below 700 hPa)
p_reg = np.linspace(p[0], 700, 10)
T_reg = np.array([])
exponent_value = (1-(7/5))/(7/5)
constant_value = (T[0] + 273.15) * (p[0] ** exponent_value)
# Calculate parcel temperature
for pressure_level in range (0, len(p_reg)):
    T_reg = np.append(T_reg, (constant_value/(p_reg[pressure_level] ** exponent_value)) - 273.15)
skew.plot(p_reg*units.hPa, T_reg*units.degC, 'ko', markerfacecolor = 'black')

#display the plot
plt.show()

# Calculating and printing Delta Z
ind = np.where(p >= 700)
index_700 = len(ind[0])
p_irr = p[0:index_700]
T_irr = np.array([])
for pLevel in range (0, len(p_irr)):
    T_irr = np.append(T_irr, (constant_value/(p_irr[pLevel] ** exponent_value)) - 273.15)
    print("At {:06.2f} hPa, Delta T: {:05.2f} C".format(p_irr[pLevel], T_irr[pLevel] - T[pLevel]))
    
# Plotting potential temperature as a function of height
thta_700 = thta[0:index_700]
height_700 = h[0:index_700]
plt.figure(1)
plt.plot(thta_700 , height_700)
plt.title('Potential Temperature vs. Height\nOUN Norman Observations at 00Z 01 July 2011')
plt.xlabel('Potential Temperature (K)')
plt.ylabel('Height (m)')
plt.show()