#Michael Johnston
#Surface Data Analysis

execfile('/home/users/xb899100/bin/projectFun.py')
import numpy as np
import matplotlib.pyplot as plt

sfcData = {'DIR' : [],
           'U' : [],
            'temperature' : [],
            'dew' : []}
dateTimes = []
with open('/home/xb899100/data/Station/SurfaceData.txt') as txtFile:
    next(txtFile)
    next(txtFile)
    for line in txtFile:
        myLine = line.split()
        DIR = int(myLine[5].split(',')[1])
        U = float(myLine[6].split(',')[0])
        temp = myLine[7].split(',')[0]
        dew = myLine[8].split(',')[0]
        
        date = myLine[4].split(',')[3]
        time = myLine[4].split(',')[4][0:2] + ':' + myLine[4].split(',')[4][2:5]
        if U < 50 and DIR < 361 and len(temp) > 0 and len(dew) > 0 and temp < '999.0' and dew < '999.0':
            sfcData['DIR'].append(DIR)
            sfcData['U'].append(U)
            sfcData['temperature'].append(float(temp))
            sfcData['dew'].append(float(dew))
            dateTimes.append(date + ' ' + time)

diurnalU = {'00': []}
diurnalDIR = {'00': []}
diurnalTemp = {'00' : []}
diurnalDew = {'00' : []}
diurnalDate = {'00' : []}
hour0 = '999'
for it in range(len(dateTimes)):
    dateTime = dateTimes[it]
    date = dateTime.split()[0]
    time = dateTime.split()[1]
    hour = time.split(':')[0]
    minutes = time.split(':')[1]
    if minutes != '00':
        hour = str(int(hour) + 1)
        if len(hour) < 2:
            hour = '0' + hour
        if hour == '24':
            hour = '00'
    
    if hour in diurnalU.keys():
        if ((minutes == '55') or (minutes == '00')) and hour0 != hour:
            diurnalU[hour].append(sfcData['U'][it])
            diurnalDIR[hour].append(sfcData['DIR'][it])
            diurnalTemp[hour].append(sfcData['temperature'][it])
            diurnalDew[hour].append(sfcData['dew'][it])
            diurnalDate[hour].append(date)
            #print('previous ' + hour0 + ' added ' + hour)
            hour0 = hour
    else:
        diurnalU[hour] = []
        diurnalDIR[hour] = []
        diurnalTemp[hour] = []
        diurnalDew[hour] = []
        diurnalDate[hour] = []
        if ((minutes == '55') or (minutes == '00')) and hour0 != hour:
            diurnalU[hour].append(sfcData['U'][it])
            diurnalDIR[hour].append(sfcData['DIR'][it])
            diurnalTemp[hour].append(sfcData['temperature'][it])
            diurnalDew[hour].append(sfcData['dew'][it])
            diurnalDate[hour].append(date)
            #print('previous ' + hour0 + ' added ' + hour)
            hour0 = hour

counts = []
for key in diurnalU.keys():
    counts.append(len(diurnalU[key]))
    
#plt.plot(diurnalU.keys(), counts, 'o')

#find the obs that are in may through october
#mean diurnal cycle
Ucycle = []
DIRcycle = []
Tcycle = []
Tdcycle = []
for hr in diurnalU.keys():
    iobs = []
    for iob in range(len(diurnalDate[hr])):
        month = diurnalDate[hr][iob][4:6]
        if month in ['05', '06', '07', '08', '09', '10']:
            iobs.append(iob)
    Ucycle.append(np.nanmean(np.array(diurnalU[hr])[iobs]))
    u, v = toComponents(np.array(diurnalU[hr])[iobs], np.array(diurnalDIR[hr])[iobs])
    umean = np.mean(u)
    vmean = np.mean(v)
    U, DIR = fromComponents(umean, vmean)
    DIRcycle.append(DIR)
    Tcycle.append(np.nanmean(np.array(diurnalTemp[hr])[iobs]))
    Tdcycle.append(np.nanmean(np.array(diurnalDew[hr])[iobs]))

#get mean wind direction for the day 1000UTC to 2300UTC
dailyU = []
dailyDIR = []
dailyDates = []
hour0 = '999'
date0 = dateTimes[0].split()[0]
timeframe = range(10, 24)
dayU = []
dayDIR = []
for it in range(len(dateTimes)):
    #only considering obs 5mins to the hour or on the hour
    dateTime = dateTimes[it]
    time = dateTime.split()[1]
    hour = time.split(':')[0]
    minutes = time.split(':')[1]
    date = dateTime.split()[0]
    #if we've moved to the next day, append and reset    
    if date != date0 and len(dayU) > 0:
        #convert to wind components to calculate mean direction
        u, v, = toComponents(np.array(dayU), np.array(dayDIR))
        umean = np.nanmean(u)
        vmean = np.nanmean(v)
        U, DIR = fromComponents(umean, vmean)
        #^gives right direction but not speed
        U = np.nanmean(dayU)
        #append
        dailyU.append(U)
        dailyDIR.append(DIR)
        dailyDates.append(datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), 23, 59, 59))
        #reset
        dayU = []
        dayDIR = []
    #round up to the next hour if not on the hour
    if minutes != '00':
        hour = str(int(hour) + 1)
        if len(hour) < 2:
            hour = '0' + hour
        if hour == '24':
            hour = '00'
    
    if int(hour) in timeframe:
        if ((minutes == '55') or (minutes == '00')) and hour != hour0:
            dayU.append(sfcData['U'][it])
            dayDIR.append(sfcData['DIR'][it])
            #print('previous ' + hour0 + ' added ' + hour)
            hour0 = hour
            date0 = date

