#read the land sea mask
#read the sounding data
#get layer mean winds
#mean by wind direction
import numpy as np
import datetime
import csv
from netCDF4 import Dataset
#land sea mask
lsm = Dataset('/glusterfs/greybls/users/xb899100/SatData/lsm-sat-coords.nc')
lsm = lsm.variables['mask'][:,:]
# land = 1, sea = 0
#change that so land = 0 and sea = 1.
lsm -= 1
lsm = np.abs(lsm)

#some common plotting variables
numberRes = int(360./myres)
sectorDIR = np.linspace(0, 360-myres, numberRes)

#array of regularly spaced heights
#Z = np.linspace(50, 15000, nlevels)
P = np.arange(1000., 99., -5.)

#read the csvs 
ndays = 2383
nlevels = len(P)
nlevelsRaw = 200
interpolated = {'Temp_p' : np.zeros((nlevels, ndays)),
          'u_p' : np.zeros((nlevels, ndays)),
            'v_p' : np.zeros((nlevels, ndays)),
            'RH_p' : np.zeros((nlevels, ndays)),
            'Z_p' : np.zeros((nlevels, ndays))}

#get the interpolated data
for key in interpolated.keys():
    with open('/home/xb899100/data/Sounding/Interpolated/' + key + '.csv', 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter = ',')
        level = 0
        for row in reader:
            interpolated[key][level,:] = row
            level += 1

#get the dates
interpolated['date']=[]
with open('/home/xb899100/data/Sounding/Interpolated/date.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        interpolated['date'].append(row[0])

#plot a histogram of the LCLp
LCLp = []
datei = 0
with open('/home/xb899100/data/Sounding/Interpolated/LCLp.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        LCLp.append(float(row[0]))
        datei += 1

MLZp = []
datei = 0
with open('/home/xb899100/data/Sounding/Interpolated/MLZp.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        MLZp.append(float(row[0]))
        datei += 1

INVTD = []
datei = 0
with open('/home/xb899100/data/Sounding/Interpolated/INVTD.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        INVTD.append(float(row[0]))
        datei += 1

CAPE = []
datei = 0
with open('/home/xb899100/data/Sounding/Interpolated/CAPE.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        CAPE.append(float(row[0]))
        datei += 1

CIN = []
datei = 0
with open('/home/xb899100/data/Sounding/Interpolated/CIN.csv', 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        CIN.append(float(row[0]))
        datei += 1

#reformat the sounding dates to match with the cloud frequency dates
soundingdate = []
for it in range(len(interpolated['date'][:])):
    year = int(float(interpolated['date'][it])*1e-06)
    month = int(float(interpolated['date'][it])*1e-04) - year*100
    day = int(float(interpolated['date'][it])*1e-02) - int(float(interpolated['date'][it])*1e-04)*100
    hour = int((float(interpolated['date'][it])*1e-02 - int(float(interpolated['date'][it])*1e-02))*100)
    
    soundingdate.append(datetime.datetime(int(year), int(month), int(day), int(hour), 00, 00))

#plot with respect to 12UTC soundings
results12 = {'Umean' : np.zeros(len(visdates)),
            'DIRmean' : np.zeros(len(visdates)),
            'visindex' : np.zeros(len(visdates))}
#sectors is the numbers corresponding to each direction
sectors = range(1, int(360/myres + 1))
#match the sounding date to the cloud frequency date and then mask it and
#calculate means for each sector and the mean sounding winds
#for the 12UTC soundings
iR = 0
iF = 0
for iS in range(len(soundingdate)):
    if iF < len(visdates) and soundingdate[iS].date() == visdates[iF].date():
        while iF < len(visdates) and soundingdate[iS].date() == visdates[iF].date():
            #iF is the satellite date that matches the sounding date in iS
            if np.min(interpolated['Z_p'][:,iS]) > 0.0:
                #checks that 'surface' altitude is higher than the highest pressure in averaging range
                U, DIR = fromComponents(interpolated['u_p'][:,iS], interpolated['v_p'][:,iS], isList = True)
                results12['Umean'][iR], results12['DIRmean'][iR] = presWeight(lowerP, upperP, p = P, U = U, DIR = DIR, units = 'm/s')
                results12['visindex'][iR] = iF
            iR += 1
            iF += 1
    elif iF < len(visdates) and soundingdate[iS].date() > visdates[iF].date():
        while iF < len(visdates) and soundingdate[iS].date() > visdates[iF].date():
            iF += 1
        while iF < len(visdates) and soundingdate[iS].date() == visdates[iF].date():
            #iF is the satellite date that matches the sounding date in iS
            if np.min(interpolated['Z_p'][:,iS]) > 0.0:
                #checks that 'surface' altitude is higher than the highest pressure in averaging range
                U, DIR = fromComponents(interpolated['u_p'][:,iS], interpolated['v_p'][:,iS], isList = True)
                results12['Umean'][iR], results12['DIRmean'][iR] = presWeight(lowerP, upperP, p = P, U = U, DIR = DIR, units = 'm/s')
                results12['visindex'][iR] = iF
            iR += 1
            iF += 1


myBins = intoBins(results12['DIRmean'], myres)
###Using station data for the bins
