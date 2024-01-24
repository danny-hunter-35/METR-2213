"""
Utilities for use in METR 3213 Atmospheric Thermodynamics
Written by Phillip Chilson
"""

import numpy as np
import urllib.request

def Gradient(z0, z1, T0, p0, rho0, lapseRate):
    # Caclulate the temperature, pressure, and density of the atmosphere based
    # a constant temperature gradient
    g = 9.81       #Acceleration of gravity (m/s/s)
    R = 287        #Gas constant for air (J/kg-K)
    
    if lapseRate == 0:
        T1 = T0
        p1 = p0*np.exp(-(g/(R*T1))*(z1-z0))      
    else:
        T1 = T0 + lapseRate*(z1 - z0)
        p1 = p0*(T1/T0)**(-g/(lapseRate*R))
    rho1 = p1/(R*T1)
    return T1, p1, rho1

def atmosphere(height):
    gamma = 1.4
    R = 287        #Gas constant for air (J/kg-K)
    
    # Initial values at critical levels (surface and layer transition)
    #Altitudes (m)
    heightInit = np.array([0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    # Pressures (Pa)
    pressInit = np.array([1.01325e5, 2.2616e+04, 5.4672e+03, 866.0263, 110.5411, 66.7019, 3.9370, 0.3711])
    # Temperatures (K)
    tempInit = np.array([288.1600, 216.6600, 216.6600, 228.6600, 270.6600, 270.6600, 214.6600, 186.9560])
    # Densities (kg/m^3)
    rhoInit = np.array([1.225, 0.3637, 0.08792, 0.013197, 0.001423, 0.0008587, 0.00006390, 0.000006917])
    #Lapse Rates for layers bounded by critical lavels (K/m)
    lapseRate = np.array([-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002])
    
    T = np.array([])
    p = np.array([])
    rho = np.array([])
    a = np.array([])
    
    for h in np.nditer(height):
        if h <= heightInit[1]:
            T1, p1, rho1 = Gradient(heightInit[0], h, tempInit[0], pressInit[0], rhoInit[0], lapseRate[0])
        elif heightInit[1] < h and h <= heightInit[2]:
            T1, p1, rho1 = Gradient(heightInit[1], h, tempInit[1], pressInit[1], rhoInit[1], lapseRate[1])
        elif heightInit[2] < h and h <= heightInit[3]:
            T1, p1, rho1 = Gradient(heightInit[2], h, tempInit[2], pressInit[2], rhoInit[2], lapseRate[2])
        elif heightInit[3] < h and h <= heightInit[4]:
            T1, p1, rho1 = Gradient(heightInit[3], h, tempInit[3], pressInit[3], rhoInit[3], lapseRate[3])
        elif heightInit[4] < h and h <= heightInit[5]:
            T1, p1, rho1 = Gradient(heightInit[4], h, tempInit[4], pressInit[4], rhoInit[4], lapseRate[4])
        elif heightInit[5] < h and h <= heightInit[6]:
            T1, p1, rho1 = Gradient(heightInit[5], h, tempInit[5], pressInit[5], rhoInit[5], lapseRate[5])
        elif heightInit[6] < h and h <= heightInit[7]:
            T1, p1, rho1 = Gradient(heightInit[6], h, tempInit[6], pressInit[6], rhoInit[6], lapseRate[6])
        else:
            T1 = np.NAN
            p1 = np.NAN
            rho1 = np.NAN
        a1 = (R*gamma*T1)**0.5
        T = np.append(T, T1)
        p = np.append(p, p1)
        rho = np.append(rho, rho1)
        a = np.append(a, a1)
    return T, p, rho, a

def getWyoming(year, month, day, hour, station):
    
    # Create the URL
    base_url = 'http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST'
    yearstr = ('&YEAR=%4.4d' % (year))
    monthstr = ('&MONTH=%2.2d' % (month))
    dayhourstr1 = ('&FROM=%2.2d%2.2d' % (day, hour))
    dayhourstr2 = ('&TO=%2.2d%2.2d' % (day, hour))
    stationstr = ('&STNM=%5.5d' % (station))
    url = base_url + yearstr + monthstr + dayhourstr1 + dayhourstr2 + stationstr
    # Create the file name
    fileName = ('%4.4d%2.2d%2.2d%2.2d_%5.5d_sounding.txt' \
                %(year, month, day, hour, station))
    # Read URL and create file
    urllib.request.urlretrieve(url, fileName)
    
    # Open the sounding file
    f = open(fileName, 'r')
    
    # Read the HTML code
    for count in range(0, 3):
        f.readline()
    # Read the sounding information
    line = f.readline()
    line = line.strip()
    line = line.replace('<H2>', '')
    line = line.replace('</H2>', '')
    soundingInfo = line
    print("%s" % (soundingInfo))
    # Read more HTML code and then the header
    for count in range(0, 5):
        f.readline()
    
    # Create the arrays
    pres = np.array([])
    hght = np.array([])
    temp = np.array([])
    dwpt = np.array([])
    relh = np.array([])
    mixr = np.array([])
    drct = np.array([])
    sknt = np.array([])
    thta = np.array([])
    thte = np.array([])
    thtv = np.array([])
        
    # Read the data
    for line in f:
        line = line.strip()
        columns = line.split()
        if columns[0] == '</PRE><H3>Station':
            break    
        pres = np.append(pres, float(columns[0]))   
        hght = np.append(hght, float(columns[1]))
        if len(columns) == 11:
            temp = np.append(temp, float(columns[2]))
            dwpt = np.append(dwpt, float(columns[3]))
            relh = np.append(relh, float(columns[4]))
            mixr = np.append(mixr, float(columns[5]))
            drct = np.append(drct, float(columns[6]))
            sknt = np.append(sknt, float(columns[7]))
            thta = np.append(thta, float(columns[8]))
            thte = np.append(thte, float(columns[9]))
            thtv = np.append(thtv, float(columns[10]))
        else:
            temp = np.append(temp, np.nan)
            dwpt = np.append(dwpt, np.nan)
            relh = np.append(relh, np.nan)
            mixr = np.append(mixr, np.nan)
            drct = np.append(drct, np.nan)
            sknt = np.append(sknt, np.nan)
            thta = np.append(thta, np.nan)
            thte = np.append(thte, np.nan)
            thtv = np.append(thtv, np.nan)
    f.close()
    return soundingInfo, pres, hght, temp, dwpt, relh, mixr, drct, sknt, thta, thte, thtv

def mean_temp_calc(T, p, p_max, p_min):
    # Approximate the mean temperature from p_max to p_min based on temperature
    # and pressure data taken from a sounding
    # Inputs
    # T: vector of temperature values (K)
    # p: vector of pressure values in any presure units
    # p_max: maximum pressure value (same units as p)
    # p_min: minimum pressure value (same units as p)
    # Ouput
    # T_mean: mean temperature
    
    # Find the indices of pressure values that fall between the target values
    # p_max and p_min
    # Find the heights of p_max-d & p_min_d and then the difference (layer
    # thickness)
    ind = np.where((p_max >= p) & (p >= p_min))[0]
    # Calculate the number of points and the actual values for p_max and p_min
    npts = len(ind)
    
    p_max = p[ind[0]]
    p_min = p[ind[-1]]
    
    # Calculate the weighted mean
    num = 0
    den = 0
    for j in np.arange(npts - 1):
        Del = np.log(p[ind[j + 1]]) - np.log(p[ind[j]])
        num = num + (T[ind[j + 1]] + T[ind[j]])/2*Del
        den = den + Del
    T_mean = num/den

    return T_mean, p_max, p_min

def hypsometric_avgT(p1, p2, T1, T2):
    # function Delta_z = hypsometric_avgT(p1, p2, T1, T2)
    # Calculate the atmospheric thickness using an approximation of the
    # hypsometric equation in which it is assumed that T_avg can be found from
    # the expression T_avg = (T1 + T2)/2
    # Inputs
    # p1 = Pressure at the first (lowest in height) level [Pa]
    # p1 = Pressure at the second (highest in height) level [Pa]
    # T1 = Temperature at the first (lowest in height) level [K]
    # T1 = Temperature at the second (highest in height) level [K]
    # Output
    # Delta_z = Layer thickness [m]
    
    Rd = 287 # J kg^-1 K^-1
    g = 9.81 # m s^-2
    
    T_avg = (T1 + T2)/2
    Delta_z = Rd*T_avg/g*np.log(p1/p2)
    return Delta_z

def get_es(T):
    e0 = 611   # reference pressure (Pa)
    T0 = 273.16   # reference temperature (K)
    Lv = 2.5e6 # latent heat of vaportization (J/kg/K)
    Rv = 461.5   # water vapor gas constant (J/kg/K)
    es = e0*np.exp(Lv/Rv*(1/T0 - 1/T))

    return es