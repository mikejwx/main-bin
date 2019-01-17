#Michael Johnston
#Deal with IGRA soundings

import numpy as np
import matplotlib.pyplot as plt
import csv

path = '/home/xb899100/OLD/Sounding/IGRA/'

#MLZp is the pressure of the top of the mixed layer found using the parcel method
#LCLp is the pressure of the lifting condensation level

igra_data = []
with open(path + 'derived/derived.txt', 'r') as igra_file:
    for line in igra_file:
        igra_data.append(line)


ndays = 2383 # number of ascents in 2012-2016 period
soundings = {'height' : np.zeros((200, ndays)), 
            'pres' : np.zeros((200, ndays)), 
            'temp' : np.zeros((200, ndays)), 
            'RH' : np.zeros((200, ndays)), 
            'uwind' : np.zeros((200, ndays)), 
            'vwind' : np.zeros((200, ndays)),
            'date' : np.zeros((1, ndays)),
            'INVTD' : np.zeros((1, ndays)),
            'MLZp' : np.zeros((1, ndays)),
            'LCLp' : np.zeros((1, ndays)),
            'CAPE' : np.zeros((1, ndays)),
            'CIN' : np.zeros((1, ndays))}
it = -1
for line in igra_data:
    myLine = line.split()
    if line[0] == '#' and int(myLine[1]) >= 2012 and int(myLine[1]) <= 2016:
        
        it += 1
        
        level = igra_data.index(line) + 1
        myInc = 0
        nlevels = int(myLine[6])
        
        soundings['date'][0,it] = int(myLine[1] + myLine[2] + myLine[3] + myLine[4])
        soundings['INVTD'][0,it] = int(line[56:61])/10.
        soundings['MLZp'][0,it] = int(line[62:67])/100.#convert to hPa
        soundings['LCLp'][0,it] = int(line[86:91])/100. #convert to hPa
        soundings['CAPE'][0,it] = int(line[146:151])
        soundings['CIN'][0,it] = int(line[152:157])
        for inc in range(nlevels):
            data_line = igra_data[level+inc]
            z = int(data_line[17:23]) #m, the calculated geopotential height
            #QC pressure
            p = int(data_line[1:7])
            #QC temperatures
            T = int(data_line[25:31])
            #QC RH
            RH = int(data_line[97:103])
            #QC wind
            u = int(data_line[113:119])
            v = int(data_line[129:135])
            
            #only include complete sounding levels
            if z > 0.:
                if p > 0.:
                    soundings['height'][myInc, it] = z
                    soundings['pres'][myInc, it] = p
                    soundings['temp'][myInc, it] = T
                    soundings['RH'][myInc, it] = RH
                    soundings['uwind'][myInc, it] = u 
                    soundings['vwind'][myInc, it] = v
                    myInc += 1

#interpolate to regular heights
execfile('/home/xb899100/bin/projectFun.py')
#Z = np.linspace(50, 15000, 300)
P = np.arange(1000., 99., -5.)*100.
zInterpolated = interpolateP(soundings['height'], P, soundings['pres']) #doesn't need to be converted to m
TempInterpolated = interpolateP(soundings['temp'], P, soundings['pres'])/10. #convert to K
uInterpolated = interpolateP(soundings['uwind'], P, soundings['pres'])/10.#convert to m/s
vInterpolated = interpolateP(soundings['vwind'], P, soundings['pres'])/10.#convert to m/s
RHInterpolated = interpolateP(soundings['RH'], P, soundings['pres'])/10. #convert to %

U_Interpolated = np.zeros_like(zInterpolated)
DIR_Interpolated = np.zeros_like(zInterpolated)
for ascent in range(it+1):
    U_Interpolated[:,ascent], DIR_Interpolated[:,ascent] = fromComponents(uInterpolated[:,ascent], vInterpolated[:,ascent], isList = True)
#save to csv
path = '/home/xb899100/data/Sounding/Interpolated/'
np.savetxt(path + 'Z_p.csv', zInterpolated, delimiter = ',')
np.savetxt(path + 'p.csv', P, delimiter = ',')
np.savetxt(path + 'Temp_p.csv', TempInterpolated, delimiter = ',')
np.savetxt(path + 'u_p.csv', uInterpolated, delimiter = ',')
np.savetxt(path + 'v_p.csv', vInterpolated, delimiter = ',')
np.savetxt(path + 'RH_p.csv', RHInterpolated, delimiter = ',')
np.savetxt(path + 'date.csv', soundings['date'][0,:], delimiter = ',')
np.savetxt(path + 'INVTD.csv', soundings['INVTD'][0,:], delimiter = ',')
np.savetxt(path + 'MLZp.csv', soundings['MLZp'][0,:], delimiter = ',')
np.savetxt(path + 'LCLp.csv', soundings['LCLp'][0,:], delimiter = ',')
np.savetxt(path + 'CAPE.csv', soundings['CAPE'][0,:], delimiter = ',')
np.savetxt(path + 'CIN.csv', soundings['CIN'][0,:], delimiter = ',')

soundingLevels = []
for day in range(ndays):
    soundingLevels.append(len(filter(lambda v: v != 0., soundings['temp'][:,day])))

