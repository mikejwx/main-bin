#!python
#Michael Johnston:
#=========================================================================
#creates a netCDF of all the cloud frequencies for each day in my dataset.
#POR 1 May to 31 October, 2012 to 2016
#=========================================================================

import numpy as np
from netCDF4 import Dataset, date2num
from os import listdir
import datetime

execfile('/home/xb899100/bin/projectFun.py')
#create a netcdf to store cloud frequency fields

dataset = Dataset('/glusterfs/greybls/users/xb899100/SatData/cloudFrequencies.nc', 'w', format = 'NETCDF4_CLASSIC')
#read a sample image
data = Dataset('/glusterfs/greybls/users/xb899100/SatData/goes13.2016.259.214519.BAND_01.nc')
#get the lat lon
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]
#get the date and time of the image
date = data.variables['imageDate'][:]

myShape = lat.shape

#create dimensions for our netCDF
x = dataset.createDimension('x', size = lon.shape[0])
y = dataset.createDimension('y', size = lon.shape[1])
time = dataset.createDimension('time', size = None)

#create variables for our netCDF
times = dataset.createVariable('time', np.float64, ('time',))
xs = dataset.createVariable('x', np.float64, ('x', 'y'))
ys = dataset.createVariable('y', np.float64, ('x', 'y'))
cldfreq = dataset.createVariable('cldfreq', np.float64, ('time', 'x', 'y'))

#set some global attributes
dataset.description = 'cloud frequencies for days in 2012-2016 averaged between 09:30:00 UTC and 23:59:59 UTC each day.'
dataset.source = 'Michael Johnston, Department of Meteorology at University of Reading, UK'

#set some variable attributes
xs.units = 'degree_east'
ys.units = 'degree_north'
cldfreq.units = 'dimensionless'
times.units = 'hours since 0001-01-01 00:00:00'
times.calendar = 'gregorian'

#writing dimension data
xs[:,:] = lon
ys[:,:] = lat

#gather the cldfreq data

thresh = 0.15

#get a list of files
files = listdir('.')
files.sort()#organize in alphanumerical order
files = files[5:-5] #ignore the first file i.e. the netCDF we're making

Lf = len(files)
data2 = []
doy0 = 122 #the first day of year 1 May 2012
visLevels = np.linspace(0.0, 1.0, 11)

timesIndex = 0
counts = []
weirdData = []
for i in range(Lf):
    #read the data
    data = Dataset('./'+files[i])
    #get the lat lon
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    myShape = lon.shape
    #get the date and time of the image
    date = data.variables['imageDate'][:]
    iTime = data.variables['imageTime'][:]

    #get the doy to check if a new day has started
    year = date/1000
    doy = date - year*1000
    
    if doy != doy0:
        counts.append(len(data2))
        #check that there isn't significant missing data
        if len(data2) > 14:
            #if we've moved onto the next day, add to the netCDF 
            
            #take the average of the cloud masks to get the frequency
            #data 2 is a list of all the cloud masks
            #define a data3 that has the frequency
            data3 = np.zeros_like(data2[0])
            
            for iData in range(len(data2)):
                if myShape == data2[iData].shape:
                    data3 += data2[iData]
                else:
                    weirdData.append(files[i])
            data3 = data3/len(data2)
            #add to the cldfreq variable in the ntCDF along the time dimension
            cldfreq[timesIndex,:,:] = data3
            
            #add to the time variable in the time dimension
            times[timesIndex] = date2num(datetime.datetime(year, 1, 1, 23, 59, 59) + datetime.timedelta(doy0 - 1), units = times.units, calendar = times.calendar)
            timesIndex += 1        
        data2 = [] #empty the data2 list, ready for the next day
        
        doy0 = doy #reset doy

    
    #calibrate the data
    #step 1 divide by 32 to convert from 16-bit to 10-bit
    excelDate = date2num(datetime.datetime(year, 1, 1, 23, 59, 59) + datetime.timedelta(doy0 - 1), units = 'days since 1900-01-01', calendar = 'gregorian')
    data = calibrateVis(data, 0.00012*excelDate - 3.72315)
    
    #step 4 convert from nominal reflectance to reflectance/albedo
    data1 = getReflectance(date, iTime, lon, lat, data)
    
    #bound the data between 0 and 1
    data1[data1 >= thresh] = 1.0
    data1[data1 < thresh] = 0.0
    
    #check that there is some sun
    mydate = str(datetime.datetime(year,1,1)+datetime.timedelta(doy - 1)).split()[0]
    hour = iTime/10000
    minute = (iTime - hour*10000)/100
    mytime = str(datetime.time(hour, minute, 00))
    SZA = np.max(getSZA(mydate, mytime, lon, lat)) #maximum solar zenith angle
    if SZA < 75:
        #add data1 to data2 list
        data2.append(data1)
    
dataset.close()


#aside: find out the data availability
year = 2012
doy = 121
counts1 = []
dates = []
count = 0
for f in range(Lf):
    if files[f].split('.')[2] == str(doy):
        count += 1
    else:
        while files[f].split('.')[2] != str(doy):
            counts1.append(count)
            dates.append(str(year)+'-'+str(doy))
            if doy > 307:
                doy = 121
                year += 1
                count = 0
            else:
                doy += 1
                count = 0
        count += 1

plt.clf()
fig = plt.gcf()
fig.set_size_inches(8, 6)
plt.plot(np.array(counts)/28.)
plt.ylabel('number of images')
plt.xlabel('date')
plt.xticks(np.arange(0, len(counts), 60), dates[0:len(counts):61], rotation = 45)
plt.title('Number of images per day, 100% = 29 images')
plt.show()