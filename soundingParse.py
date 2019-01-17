#read the land sea mask
#read the sounding data
#get layer mean winds
#mean by wind direction
import numpy as np
import datetime
import csv
from netCDF4 import Dataset, num2date, date2num

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
#Z = np.arange(50, 15001, 50.)
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

#get the sounding index that matches the satellite image dates
dateindexes00 = []
dateindexes12 = []
for date in clddates:
    for i, j in enumerate(soundingdate):
        if (j - date) == datetime.timedelta(-1, 43201):
            dateindexes12.append(i)
        elif (j - date) == datetime.timedelta(-1, 1):
            dateindexes00.append(i)

#can now eitherplot with respect to 12UTC or 00UTC soundings
results00 = {'Umean' : np.zeros(len(dateindexes00)),
            'DIRmean' : np.zeros(len(dateindexes00)),
            'cldindex' : np.zeros(len(dateindexes00))}
results12 = {'Umean' : np.zeros(len(dateindexes12)),
            'DIRmean' : np.zeros(len(dateindexes12)),
            'cldindex' : np.zeros(len(dateindexes12))}
#sectors is the numbers corresponding to each direction
sectors = range(1, int(360/myres + 1))
#match the sounding date to the cloud frequency date and then mask it and
#calculate means for each sector and the mean sounding winds
#for the 00UTC soundings
iR = 0
for iS in dateindexes00:
    #iS is the sounding index that has a matching satellite image
    iF = np.min(np.nonzero((soundingdate[iS] - clddates) == datetime.timedelta(-1, 1)))
    #iF is the satellite date that matches the sounding date in iS
    if P[0] > lowerP:
        if np.min(interpolated['Z_p'][:,iS]) > 0.0:
        #checks that surface altitude is higher than the highest pressure in averaging range
            U, DIR = fromComponents(interpolated['u_p'][:,iS], interpolated['v_p'][:,iS], isList = True)
            results00['Umean'][iR], results00['DIRmean'][iR] = presWeight(lowerP, upperP, p = P, U = U, DIR = DIR, units = 'm/s')
            results00['cldindex'][iR] = iF
            iR += 1

#for the 12UTC soundings
iR = 0
for iS in dateindexes12:
    #iS is the sounding index that has a matching satellite image
    iF = np.min(np.nonzero((soundingdate[iS] - clddates) == datetime.timedelta(-1, 43201)))

    #iF is the satellite date that matches the sounding date in iS
    if P[0] > lowerP:
        if np.min(interpolated['Z_p'][:,iS]) > 0.0:
            #checks that 'surface' altitude is higher than the highest pressure in averaging range
            U, DIR = fromComponents(interpolated['u_p'][:,iS], interpolated['v_p'][:,iS], isList = True)
            results12['Umean'][iR], results12['DIRmean'][iR] = presWeight(lowerP, upperP, p = P, U = U, DIR = DIR, units = 'm/s')
            results12['cldindex'][iR] = iF
            iR += 1

myBins = intoBins(results12['DIRmean'], myres)