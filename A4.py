# Load the packages
import numpy as np
import matplotlib.pyplot as plt
# ---------------------------------------------------------------------
# Constant parameters
# ---------------------------------------------------------------------
# Parameters for dry air
R_d = 287      # Gas constant (J/kg/K)
c_p = 1005     # Specific heat capacity (J/kg/K)
k_d = R_d/c_p  # Dimensionless
# Parameters for water
R_v = 461      # Gas constant (J/kg/K)
L_v = 2.5e6    # Latent heat of vaporization (J/kg)
e_0 = 611       # Triple point pressure (Pa)
T_0 = 273       # Triple point temperature (K)
# Other
epsilon = 0.622

# ---------------------------------------------------------------------
# Parameters used in the example mentioned above
# ---------------------------------------------------------------------
# Pressure is in Pa and temerature is in K
p_beg = 100000
T_beg = 293
Td_beg = 288
p_step = 100
p_top = 70000
p_end = 95000

# Assign the pressures for parcel during ascent
pup_p = np.arange(p_beg, p_top - p_step, -p_step)

# Calculate the vapor pressure and mixing ratio based on the surface conditions, e.g.,
e = e_0*np.exp(L_v/R_v*(1/T_0 - 1/Td_beg))
w = epsilon*e/(p_beg - e)

# ---------------------------------------------------------------------
# Calculate the temperature and dew point temperature for the ascent
# from the surface up to the LCL (or just below it) using something like
# ---------------------------------------------------------------------
# Set up the arrays for temperature and dew-point temperature
pup_T = np.zeros(pup_p.size)
pup_Td = np.zeros(pup_p.size)

# surfact to LCL
for i in range(pup_p.size):
    T = T_beg*(pup_p[i]/p_beg)**k_d
    e = w/(w + epsilon)*pup_p[i]
    Td = (1/T_0 - R_v/L_v*np.log(e/e_0))**(-1)
    if T < Td:
        ind = i - 1
        break
    else:
        pup_T[i] = T
        pup_Td[i] = Td

# Initial pressure and temperature of LCL
T_LCL = pup_T[ind]
print(T_LCL)
p_LCL = pup_p[ind]
# ---------------------------------------------------------------------
# Calculate the temperature and dew point temperature for the ascent
# from the LCL (or just above it) to the top using something like
# ---------------------------------------------------------------------
for i in range(ind + 1, pup_p.size):
    e_s = e_0*np.exp(L_v/R_v*(1/T_0 - 1/pup_T[i-1]))
    w_s = epsilon*(e_s/(pup_p[i-1] - e_s))
    dTdp = ((R_d/c_p)*pup_T[i-1] + (L_v/c_p)*w_s)/(pup_p[i-1]*(1 + (w_s * epsilon * L_v**2)/(c_p * R_d * T**2)))
    pup_T[i] = pup_T[i-1] - dTdp*p_step
    pup_Td[i] = pup_T[i]
    
# ---------------------------------------------------------------------
# Calculate the temperature and dew point temperature for the descent
# from the top to the ending pressure level using something like
# ---------------------------------------------------------------------
# assign the value of the temperature at the top from the last value in pup_T
T_top = pup_T[-1]
# Create arrays for pressure and temperature during descent
pdn_p = np.arange(p_top, p_end + p_step, p_step)
pdn_T = np.zeros(pdn_p.size)
# calculate the temperature on the descent from the Poisson's equation
for i in range(pdn_p.size):
    pdn_T[i] = T_top*(pdn_p[i]/p_top)**k_d

# ---------------------------------------------------------------------
# Plot the results using something like
# ---------------------------------------------------------------------
# fig = plt.figure(figsize=(9, 9))
# plt.plot(pup_T - 273, pup_p/100, '-r')
# plt.plot(pdn_T - 273, pdn_p/100, '-b')
# plt.gca().invert_yaxis()
# add title and axes labels
