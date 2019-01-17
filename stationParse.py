#reads the surface winds measured on each runway at 10minute resolution.
#results in myBins to be useable with CTID_vis.py
execfile('/home/users/xb899100/bin/projectFun.py')
data12 = {'spd' : [],
        'DIR' : [],
        'date-time' : []}

for year in range(2012, 2017):
    with open('/home/xb899100/data/Station/WIND_KT_SITE_12_'+str(year)+'.txt') as txtFile:
        next(txtFile)
        for line in txtFile:
            myLine = line.split()
            if len(myLine) == 9:
                date = myLine[0]
                time = myLine[1]
                spd = myLine[6]
                DIR = myLine[7]
                data12['spd'].append(float(spd)*0.5144)
                data12['DIR'].append(int(DIR))
                data12['date-time'].append(date +' '+time)

#data30 = {'spd' : [],
#        'DIR' : [],
#        'date-time' : []}
#
#for year in range(2012, 2017):
#    with open('/home/xb899100/data/Station/WIND_KT_SITE_30_'+str(year)+'.txt') as txtFile:
#        next(txtFile)
#        for line in txtFile:
#            myLine = line.split()
#            if len(myLine) == 9:
#                date = myLine[0]
#                time = myLine[1]
#                spd = myLine[6]
#                DIR = myLine[7]
#                data30['spd'].append(float(spd)*0.5144)
#                data30['DIR'].append(int(DIR))
#                data30['date-time'].append(date +' '+time)

myDIRs = data12['DIR']
mySpds = data12['spd']
myDates = data12['date-time']

#months = ['05', '06', '07', '08', '09', '10']
#hours = ['09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']
#minutes = ['20', '50']
#for date in data12['date-time']:
#    month = date.split('-')[1]
#    hour = date.split()[1].split(':')[0]
#    minute = date.split()[1].split(':')[1]
#    if (month in months) and (hour in hours) and (minute in minutes):
#        my_i = np.where(np.array(data12['date-time']) == date)[0]
#        myDIRs.append(data12['DIR'][my_i])
#        mySpds.append(data12['spd'][my_i])
#        myDates.append(data12['date-time'][my_i])
if 'myres' in globals():
    myBins = intoBins(np.array(myDIRs), myres)
else:
    print "WARNING: 'myBins' not created! 'myres' is not defined!"