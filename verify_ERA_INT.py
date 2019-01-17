#Code to compare soundings to ERA-interim data
#focusing on temperature, specific humidity, U, and V-winds

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2num
import datetime
import math
import csv
execfile('/home/xb899100/bin/projectFun.py')
execfile('/home/xb899100/bin/SkewT.py')
execfile('/home/xb899100/main/ERA-int/Read_ERA.py')

def get_highs(x, y, z, window = 5):
    #finds the local maxima
    #returns the lon, lat, and value of all the local maxima
    #uses a default five point window
    window = int(window)
    if window%2 == 0:
        window += 1
    w0 = window/2
    w1 = window - w0
    x_max = []
    y_max = []
    z_max = []
    for i in np.arange(w0, len(x) - w1):
        for j in np.arange(w0, len(y) - w1):
            if z[j,i] == np.max(z[(j-w0):(j+w1),(i-w0):(i+w1)]):
                x_max.append(x[i])
                y_max.append(y[j])
                z_max.append(z[j,i])
    return x_max, y_max, z_max

def get_lows(x, y, z, window = 5):
    #finds the local maxima
    #returns the lon, lat, and value of all the local maxima
    #uses a default five point window
    window = int(window)
    if window%2 == 0:
        window += 1
    w0 = window/2
    w1 = window - w0
    x_min = []
    y_min = []
    z_min = []
    for i in np.arange(w0, len(x) - w1):
        for j in np.arange(w0, len(y) - w1):
            if z[j,i] == np.min(z[(j-w0):(j+w1),(i-w0):(i+w1)]):
                x_min.append(x[i])
                y_min.append(y[j])
                z_min.append(z[j,i])
    return x_min, y_min, z_min

plt_next = True
time = 0
while plt_next:
    plt.clf()
    plt.contour(lon-360, lat, mslp[time,:,:], cmap = 'Blues',levels = np.arange(950, 1050, 2)*100.)
    plt.grid()
    x0, y0, z0 = get_highs(lon, lat, mslp[time,:,:], 3)
    for H in range(len(x0)):
        plt.plot(x0[H]-360., y0[H], 'bo')
        plt.text(x0[H]-360., y0[H], 'H\n' + str(round(z0[H]/100.)) + 'hPa', color = 'b')
    x1, y1, z1 = get_lows(lon, lat, mslp[time,:,:], 3)
    for L in range(len(x1)):
        plt.plot(x1[L]-360., y1[L], 'ro')
        plt.text(x1[L]-360., y1[L], 'L\n' + str(round(z1[L]/100.)) + 'hPa', color = 'r')
    plt.title(str(num2date(times[time], times_units, times_calendar)))
    plt.show()
    time += 1
    plt_next = int(raw_input('Next?') or True)
    
#Bermuda is near 32.25 and 295.5 in the data
BDA_x = np.where(lon == 295.5)[0][0]
BDA_y = np.where(lat == 32.25)[0][0]

#read in the sounding data
execfile('/home/xb899100/bin/readSoundingData.py')
#most data in the variable 'interpoalted'

for time_index in range(len(times)):
    #convert ERA-Int time to the same datetime.datetime format as the soudningdate
    time_convert = num2date(times[time_index], times_units, times_calendar)
    sounding_index = np.where(np.array(soundingdate) == time_convert)[0]
    if sounding_index > 0:
        #plot the ERA skew-T
        t0 = 273.15
        ERA_T = temperatures[time_index,:,BDA_y, BDA_x]-t0
        ERA_T = ERA_T[::-1]
        ERA_P = p_level[::-1]
        ERA_Td = getDew(hum[time_index,:,BDA_y,BDA_x], p_level)-t0
        ERA_Td = ERA_Td[::-1]
        ERA_u = u[time_index, :, BDA_y, BDA_x]/0.5144
        ERA_u = ERA_u[::-1]
        ERA_v = v[time_index, :, BDA_y, BDA_x]/0.5144
        ERA_v = ERA_v[::-1]
        #plot the OBS skew-T
        OBS_T = interpolated['Temp_p'][:,int(sounding_index)]-t0
        OBS_Td = getDew(getQ(OBS_T,interpolated['RH_p'][:,int(sounding_index)],P),P)-t0
        OBS_u = interpolated['u_p'][:, int(sounding_index)]/0.5144
        OBS_v = interpolated['v_p'][:, int(sounding_index)]/0.5144
        #plot both on the same Skew-T
        plt.clf()
        gs = gridspec.GridSpec(1, 2, width_ratios = [4, 1])
        gs.update(wspace=0.0, hspace=0.0)
        fig = plt.gcf()
        fig.set_size_inches(6, 8)
        paper()
        plt.plot(skew(OBS_T, P), P, color = 'r', label = 'OBS Temp')
        plt.plot(skew(OBS_Td, P), P, color = 'g', label = 'OBS Dew')
        #plt.plot(skew(ERA_T, ERA_P), ERA_P, 'r--', label = 'ERA Temp')
        #plt.plot(skew(ERA_Td, ERA_P), ERA_P, 'g--', label = 'ERA Dew')
        plt.legend(loc = 2)
        plt.title(str(time_convert)+' UTC')
        col = 'cyan'
        alpha1 = 1.0
        #windBarbs(ERA_u, ERA_v, p_level)
        #plt.text(0, 70, 'ERA wind', color = col)
        col = 'k'
        alpha1 = 1.0
        windBarbs(OBS_u, OBS_v, P)
       #plt.text(-0.5, 90, 'OBS wind', color = col, alpha = alpha1)
        plt.show()
        myNext = raw_input('Next sounding?') or 'y'
        if myNext != 'y':
            break

#calculate some summary statistics on the errors taking obs to be the truth
n_soundings = 0
t0 = 273.15
for time_index in range(len(times)):
    #convert ERA-Int time to the same datetime.datetime format as the soudningdate
    time_convert = num2date(times[time_index], times_units, times_calendar)
    sounding_index = np.where(np.array(soundingdate) == time_convert)[0]
    if sounding_index > 0:
        #OBS skew-T
        OBS_T1 = interpolated['Temp_p'][:,int(sounding_index)]-t0
        OBS_Td1 = getDew(getQ(OBS_T1,interpolated['RH_p'][:,int(sounding_index)],P),P)-t0
        OBS_u1 = interpolated['u_p'][:, int(sounding_index)]/0.5144
        OBS_v1 = interpolated['v_p'][:, int(sounding_index)]/0.5144
        if min(OBS_T1) != -t0:
            #ERA skew-T
            ERA_T1 = temperatures[time_index,:,BDA_y, BDA_x]-t0
            ERA_T1 = ERA_T1[::-1]
            ERA_P = p_level[::-1]
            ERA_Td1 = getDew(hum[time_index,:,BDA_y,BDA_x], p_level)-t0
            ERA_Td1 = ERA_Td1[::-1]
            ERA_u1 = u[time_index, :, BDA_y, BDA_x]/0.5144
            ERA_u1 = ERA_u1[::-1]
            ERA_v1 = v[time_index, :, BDA_y, BDA_x]/0.5144
            ERA_v1 = ERA_v1[::-1]
            
            #get the OBS and ERA onto same grid
            #neglect all but the common p_levels
            ERA_T = []
            ERA_Td = []
            ERA_u = []
            ERA_v = []
            OBS_T = []
            OBS_Td = []
            OBS_u = []
            OBS_v = []
            new_pressure = []
            for pressure in ERA_P:
                if pressure >= 100.:
                    ERA_index = int(np.where(ERA_P == pressure)[0])
                    OBS_index = int(np.where(P == pressure)[0])
                    if pressure in P:
                        new_pressure.append(pressure)
                        ERA_T.append(ERA_T1[ERA_index])
                        ERA_Td.append(ERA_Td1[ERA_index])
                        ERA_u.append(ERA_u1[ERA_index])
                        ERA_v.append(ERA_v1[ERA_index])
                        OBS_T.append(OBS_T1[OBS_index])
                        OBS_Td.append(OBS_Td1[OBS_index])
                        OBS_u.append(OBS_u1[OBS_index])
                        OBS_v.append(OBS_v1[OBS_index])
            ERA_T = np.array(ERA_T)
            ERA_Td = np.array(ERA_Td)
            ERA_u = np.array(ERA_u)
            ERA_v = np.array(ERA_v)
            OBS_T = np.array(OBS_T)
            OBS_Td = np.array(OBS_Td)
            OBS_u = np.array(OBS_u)
            OBS_v = np.array(OBS_v)
            if n_soundings == 0:
                temp_error = ERA_T - OBS_T
                dew_error = ERA_Td - OBS_Td
                u_error = ERA_u - OBS_u
                v_error = ERA_v - OBS_v
            else:
                temp_error = np.vstack((temp_error, ERA_T - OBS_T))
                dew_error = np.vstack((dew_error, ERA_Td - OBS_Td))
                u_error = np.vstack((u_error, OBS_u - ERA_u))
                v_error = np.vstack((v_error, OBS_v - ERA_v))
            n_soundings += 1

mean_temp_error = []
mean_dew_error = []
mean_u_error = []
mean_v_error = []
for ilevel in range(len(temp_error[0,:])):
    mean_temp_error.append(np.nanmean(temp_error[:,ilevel]))
    mean_dew_error.append(np.nanmean(dew_error[:,ilevel]))
    mean_u_error.append(np.nanmean(u_error[:,ilevel]))
    mean_v_error.append(np.nanmean(v_error[:,ilevel]))

fig, ax = plt.subplots()
fig.set_size_inches(20, 5)
plt.subplot(141)
plt.semilogy(mean_temp_error, new_pressure)
plt.gca().invert_yaxis()
plt.yticks(wmo_p, wmo_p)
plt.grid()
plt.title('Mean Temperature Error')
plt.xlabel('deg C')

plt.subplot(142)
plt.semilogy(mean_dew_error, new_pressure)
plt.gca().invert_yaxis()
plt.yticks(wmo_p, wmo_p)
plt.grid()
plt.title('Mean Dew Point Error')
plt.xlabel('deg C')

plt.subplot(143)
plt.semilogy(mean_u_error, new_pressure)
plt.gca().invert_yaxis()
plt.yticks(wmo_p, wmo_p)
plt.grid()
plt.title('Mean U-vector Error')
plt.xlabel('kts')

plt.subplot(144)
plt.semilogy(mean_v_error, new_pressure)
plt.gca().invert_yaxis()
plt.yticks(wmo_p, wmo_p)
plt.grid()
plt.title('Mean V-vector Error')
plt.xlabel('kts')