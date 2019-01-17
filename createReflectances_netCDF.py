#!python
#Michael Johnston:
###Code to read in netCDFs of GOES 13 visible satellite data, then save all the
###individual images into one netCDF

import numpy as np
from netCDF4 import Dataset, date2num
from os import listdir
import datetime

execfile('/home/xb899100/bin/projectFun.py')
read_path = '/glusterfs/greybls/users/xb899100/SatData/'
save_path = '/home/xb899100/data/'
save_file = 'testData.nc' #'reflectances.nc'
#create a netcdf to store reflectance fields
dataset    = Dataset(save_path + save_file, 'w', format = 'NETCDF4_CLASSIC')
#read a sample image
data       = Dataset(read_path + 'goes13.2016.259.214519.BAND_01.nc')
#get the lat lon
lat        = data.variables['lat'][:]
lon        = data.variables['lon'][:]
#get the date and time of the image
date       = data.variables['imageDate'][:]

myShape    = lat.shape

#create dimensions for our netCDF
x          = dataset.createDimension('x', size = lon.shape[0])
y          = dataset.createDimension('y', size = lon.shape[1])
time       = dataset.createDimension('time', size = None)

#create variables for our netCDF
times      = dataset.createVariable('time', np.float64, ('time',))
xs         = dataset.createVariable('x', np.float64, ('x', 'y'))
ys         = dataset.createVariable('y', np.float64, ('x', 'y'))
albedo     = dataset.createVariable('albedo', np.float64, ('time', 'x', 'y'))

#set some global attributes
dataset.description  = 'cloud frequencies for days in May through October 2012-2016.'
dataset.source       = 'Michael Johnston, Department of Meteorology at University of Reading, UK'

#set some variable attributes
xs.units       = 'degree_east'
ys.units       = 'degree_north'
albedo.units   = 'dimensionless'
times.units    = 'hours since 0001-01-01 00:00:00'
times.calendar = 'gregorian'

#writing dimension data
xs[:,:] = lon
ys[:,:] = lat

#get a list of files
files = listdir('.')
files.sort()        #organize in alphanumerical order
files = files[7:-6] #ignore the definitely not satellite image files

Lf        = len(files)
data2     = []
doy0      = 122     #the first day of year 1 May 2012

timesIndex = 0
for i in range(144):
    test = files[i].split('.')
    #check that it is a netCDF file and that it is a goes image
    if len(test) == 6 and (test[-1] == 'nc'):       
        #read the data
        data = Dataset('./'+files[i])
        #get the lat lon
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        #get the date and time of the image
        date = data.variables['imageDate'][:]
        iTime = data.variables['imageTime'][:]
    
        #get the doy to check if a new day has started
        year = date/1000
        doy = date - year*1000
        
        #calibrate the data
        #steps 1-3 divide by 32 to convert from 16-bit to 10-bit, calibrate and post-launch correct
        excelDate = date2num(datetime.datetime(year, 1, 1, 23, 59, 59) + datetime.timedelta(doy), units = 'days since 1900-01-01', calendar = 'gregorian')
        data = calibrateVis(data, 0.00012*excelDate - 3.72315)
        
        #step 4 convert from nominal reflectance to reflectance/albedo
        data1 = getReflectance(date, iTime, lon, lat, data)
        
        #add to the albedo variable in the netCDF along the time dimension
        albedo[timesIndex,:,:] = data1[0:myShape[0], 0:myShape[1]]
        
        #add to the time variable in the time dimension
        times[timesIndex] = date2num(datetime.datetime(year, 1, 1, iTime/10000, iTime/100 - 100*(iTime/10000), 0) + datetime.timedelta(doy-1), units = times.units, calendar = times.calendar)
        timesIndex += 1
    
dataset.close()


#aside: find out the data availability
#year = 2012
#doy = 121
#counts = []
#dates = []
#count = 0
#for f in range(Lf):
#    if files[f].split('.')[2] == str(doy):
#        count += 1
#    else:
#        while files[f].split('.')[2] != str(doy):
#            counts.append(count)
#            dates.append(str(year)+'-'+str(doy))
#            if doy > 307:
#                doy = 121
#                year += 1
#                count = 0
#            else:
#                doy += 1
#                count = 0
#        count += 1

#plt.clf()
#fig = plt.gcf()
#fig.set_size_inches(8, 6)
#plt.plot(counts)
#plt.ylabel('number of images')
#plt.xlabel('date')
#plt.xticks(np.arange(0, len(counts), 60), dates[0:len(counts):61], rotation = 45)
#plt.title('Number of images per day, 100% = 29 images')