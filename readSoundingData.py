#read sounding data
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date, date2num
import datetime
import math
import csv

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
