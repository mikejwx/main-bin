#Michael Johnston
#Surface Data Analysis

execfile('/home/xb899100/bin/projectFun.py')

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
        hour = myLine[4].split(',')[4][0:2]
        minute = myLine[4].split(',')[4][2:5]
        
        if hour in ['08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']:
            if minute == '00':
                time = hour + ':' + minute
                if U < 50 and DIR < 361 and len(temp) > 0 and len(dew) > 0 and temp < '999.0' and dew < '999.0':
                    sfcData['DIR'].append(DIR)
                    sfcData['U'].append(U)
                    sfcData['temperature'].append(float(temp))
                    sfcData['dew'].append(float(dew))
                    dateTimes.append(date + ' ' + time)
            elif minute == '55':
                hour = str(int(hour) + 1)
                if len(hour) < 2:
                    hour = '0'+hour
                time = hour + ':' + '00'
                if U < 50 and DIR < 361 and len(temp) > 0 and len(dew) > 0 and temp < '999.0' and dew < '999.0':
                    sfcData['DIR'].append(DIR)
                    sfcData['U'].append(U)
                    sfcData['temperature'].append(float(temp))
                    sfcData['dew'].append(float(dew))
                    dateTimes.append(date + ' ' + time)

myBins1 = intoBins(np.array(sfcData['DIR']), myres)
myBins = []
for mydate in visdates:
    #puts the date of the satellite image into the format of the surface data
    date = str(mydate.date()).split('-')
    date = date[0] + date[1] + date[2]
    #round the time of the satellite image to the nearest hour
    time = str(mydate.time()).split(':')
    if time[1] == '15':
        time = time[0] + ':00'
    else:
        hour = int(time[0]) + 1
        if len(str(hour)) < 2:
            hour = '0' + str(hour)
        time = str(hour) + ':00'
    mydateTime = date + ' ' + time
    if mydateTime in dateTimes:
        myBins.append(myBins1[dateTimes.index(mydateTime)])
    else:
        myBins.append(0.0)