#All the IGRA data for my da test
import numpy as np
import matplotlib.pyplot as plt
import csv
import datetime as dt

path = '/home/xb899100/OLD/Sounding/IGRA/'

#MLZp is the pressure of the top of the mixed layer found using the parcel method
#LCLp is the pressure of the lifting condensation level

igra_data = []
with open(path + 'derived/derived.txt', 'r') as igra_file:
    for line in igra_file:
        igra_data.append(line)

#a dictionary storing lists of arrays of variables
soundings = {'pres' : [], 
            'temp' : [], 
            'RH' : [], 
            'uwind' : [], 
            'vwind' : [],
            'date' : []}
for line in igra_data:
    myLine = line.split()
    if line[0] == '#':
        level = igra_data.index(line) + 1
        nlevels = int(myLine[6])
        
        #store the date of this sounding
        soundings['date'].append(dt.datetime(int(myLine[1]), int(myLine[2]), int(myLine[3]), int(myLine[4]), 0, 0))
        #initialise arrays for the variables of interest
        pres = np.array([])
        temp = np.array([])
        RH = np.array([])
        uwind = np.array([])
        vwind = np.array([])
        for inc in range(nlevels):
            data_line = igra_data[level+inc]
            z = int(data_line[17:23]) #m, the calculated geopotential height
            #QC pressure
            p = int(data_line[1:7])
            #QC temperatures
            T = int(data_line[25:31])
            #QC RH
            rh = int(data_line[97:103])
            #QC wind
            u = int(data_line[113:119])
            v = int(data_line[129:135])
            
            #only include complete sounding levels
            if z > 0.:
                if p > 0.:
                    #concatenate the new obs onto those variables' arrays
                    pres = np.concatenate((pres, np.array([p])), 0)
                    temp = np.concatenate((temp, np.array([T])), 0)
                    RH = np.concatenate((RH, np.array([rh])), 0)
                    uwind = np.concatenate((uwind, np.array([u])), 0)
                    vwind = np.concatenate((vwind, np.array([v])), 0)
        #append the new arrays the lists of variables
        soundings['pres'].append(pres)
        soundings['temp'].append(temp)
        soundings['RH'].append(RH)
        soundings['uwind'].append(uwind)
        soundings['vwind'].append(vwind)

### Some summary plots
#my_time = 12
#sfc_temps = [soundings['temp'][ob][0] for ob in range(len(soundings['date'])) if soundings['date'][ob].time() == dt.datetime(1900, 1, 1, my_time).time()]
#sfc_dates = [soundings['date'][ob] for ob in range(len(soundings['date'])) if soundings['date'][ob].time() == dt.datetime(1900, 1, 1, my_time).time()]
#
#plt.plot(sfc_dates, np.array(sfc_temps)/10.-273.15, 'ko')
#plt.plot(sfc_dates[sfc_temps.index(np.min(sfc_temps))], np.min(sfc_temps)/10 - 273.15, 'bo')
#plt.plot(sfc_dates[sfc_temps.index(np.max(sfc_temps))], np.max(sfc_temps)/10 - 273.15, 'ro')
#plt.title(str(my_time) + 'UTC min: '+str(np.min(sfc_temps)/10.-273.15) + ', ' + str(sfc_dates[sfc_temps.index(np.min(sfc_temps))].date()) + 
#'\n' + str(my_time) + 'UTC max: ' + str(np.max(sfc_temps)/10.-273.15) + ', ' + str(sfc_dates[sfc_temps.index(np.max(sfc_temps))].date()))
#plt.show()
#
#sfc_windspd = [np.sqrt(soundings['uwind'][ob][0]**2 + soundings['vwind'][ob][0]**2) for ob in range(len(soundings['date']))]
#all_sfc_dates = [soundings['date'][ob] for ob in range(len(soundings['date']))]
#sfc_windspd = [ob if ob < 500. else np.nan for ob in sfc_windspd]
#
#plt.plot(all_sfc_dates, np.array(sfc_windspd)/0.5144/10., 'ko')
#sfc_pres = [soundings['pres'][ob][0] if soundings['pres'][ob][0] != 100000. else soundings['pres'][ob][1] for ob in range(len(soundings['date']))]
#
#plt.plot(all_sfc_dates, np.array(sfc_pres)/100., 'ko')