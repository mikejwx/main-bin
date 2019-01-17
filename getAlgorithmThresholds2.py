#code to calculate the optimised alpha and beta for the algorithm based on a manual classification.
execfile('/home/xb899100/bin/projectFun.py') #main project functions
execfile('/glusterfs/greybls/users/xb899100/SatData/readReflectances.py')
import csv
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt
from scipy import stats

print 'getAlgorithmThresholds2 -> Reading Manual Classifications'
#read in the manual classification
manual_CT = []
manual_NT = []
manual_OB = []

###old manual classification is in training_XX.csv, new manual classification is in manual_XX.csv
#cloud trail dates
with open('/home/users/xb899100/main/CloudTrailVis_Testing/Manual/Manual_CT.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter = ',')
    for row in spamreader:
        manual_CT.append(row)
#non-trail dates
with open('/home/users/xb899100/main/CloudTrailVis_Testing/Manual/Manual_NT.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter = ',')
    for row in spamreader:
        manual_NT.append(row)
#obscured dates
with open('/home/users/xb899100/main/CloudTrailVis_Testing/Manual/Manual_OB.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter = ',')
    for row in spamreader:
        manual_OB.append(row)

#get rid of weird indexing
manual_CT = manual_CT[0]
manual_NT = manual_NT[0]
manual_OB = manual_OB[0]

if getAlphaMethod[experiment] == 1:
    mask = Mask2(radius = outerR, lon = lon, lat = lat, land = lsmQ, inner = innerR, res = myres) # <-circular mask
elif getAlphaMethod[experiment] == 2: 
    mask = Mask3(radius = outerR, lon = lon, lat = lat, res = myres, inner_r = innerR) # <- buffer mask

print 'getAlgorithmThresholds2 -> Finding the optimal value for alpha'
#Alpha - the cloud fraction in the radius around Bermuda.
cloud_fractions_NO = []
cloud_fractions_OB = []

#Do for non-obscured scenes
for date in manual_CT+manual_NT:
    if date != '':
        #get the image
        dt_date = dt.strptime(date, '%Y-%m-%d %H:%M:%S')
        image = getOneFrequency(dt_date)
        
        #convert to a cloud mask
        cloud_mask = image*1.0
        cloud_mask[image > cldmask] = 1.0
        cloud_mask[image <= cldmask] = 0.0
        
        cloud_fractions_NO.append(np.nanmean(cloud_mask*mask/mask))

#do for obscured scenes
for date in manual_OB:
    if date != '':
        #get the image
        dt_date = dt.strptime(date, '%Y-%m-%d %H:%M:%S')
        image = getOneFrequency(dt_date)
        
        #convert to a cloud mask
        cloud_mask = image*1.0
        cloud_mask[image > cldmask] = 1.0
        cloud_mask[image <= cldmask] = 0.0
        
        cloud_fractions_OB.append(np.nanmean(cloud_mask*mask/mask))

NO = plt.hist(cloud_fractions_NO, bins = np.arange(0, 1.01, 0.01), cumulative = -1, normed = True, label = 'Non-OB dist')
OB = plt.hist(cloud_fractions_OB, bins = np.arange(0, 1.01, 0.01), cumulative = 1, normed = True, alpha = 0.65, label = 'OB dist')
alpha = NO[1][np.where(np.abs(NO[0] - OB[0]) == np.nanmin(np.abs(NO[0] - OB[0])))[0][0]]
plt.show()

print 'getAlgorithmThresholds2 -> alpha = ' + str(alpha)

#define methods for calculating the upwind and downwind cloud fraction
def get_upwind(method, cm, myBin):
    "methods for getting the upwind value"
    "method 1 = use the exact upwind sector"
    "method 2 = use the maximum of the upwind sectors in a window"
    "method 3 = use the mean cloud fraction in the window"
    if method == 1:
        upwind = np.nanmean(cm[mask == myBin])
    elif method == 2:
        upwinds = []
        adjustment = int(window/myres)
        for adj in range(-adjustment, 1+adjustment):
            upBin = myBin + adj
            #prevent direction folding
            if upBin < 1:
                upBin += numberRes
            elif upBin > numberRes:
                upBin -= numberRes
            #get the cloud fraction in the upwind direction
            upwind = np.nanmean(cm[mask == upBin])
            upwinds.append(upwind)
        upwind = np.nanmax(upwinds)
    elif method == 3:
        upwinds = []
        adjustment = int(window/myres)
        for adj in range(-adjustment, 1+adjustment):
            upBin = myBin + adj
            #prevent direction folding
            if upBin < 1:
                upBin += numberRes
            elif upBin > numberRes:
                upBin -= numberRes
            #get the cloud fraction in the upwind window
            upwinds = upwinds + list(cm[mask == upBin])
        upwind = np.nanmean(np.array(upwinds))
    return upwind

def get_downwind(method, cm, myBin):
    "methods for getting the downwind value"
    "method 1 = use the exact downwind sector"
    "method 2 = use the maximum of the downwind sectors in a window"
    "method 3 = use the mean cloud fraction in the window"
    if method == 1:
        downBin = myBin - numberRes/2
        if downBin < 1:
            #prevent direction folding
            downBin += numberRes
        downwind = np.nanmean(cm[mask == downBin])
    elif method == 2:
        downwinds = []
        adjustment = int(window/myres)
        for adj in range(-adjustment, 1+adjustment):
            downBin = myBin - numberRes/2 + adj
            #prevent direction folding
            if downBin < 1:
                downBin += numberRes
            elif downBin > numberRes:
                downBin -= numberRes
            #get the cloud fraction in the upwind direction
            downwind = np.nanmean(cm[mask == downBin])
            downwinds.append(downwind)
        downwind = np.nanmax(downwinds)
    elif method == 3:
        downwinds = []
        adjustment = int(window/myres)
        for adj in range(-adjustment, 1+adjustment):
            downBin = myBin - numberRes/2 + adj
            #prevent direction folding
            if downBin < 1:
                downBin += numberRes
            elif downBin > numberRes:
                downBin -= numberRes
            #get the cloud fraction in the upwind window
            downwinds = downwinds + list(cm[mask == downBin])
        downwind = np.nanmean(np.array(downwinds))
    return downwind

print 'getAlgorithmThresholds2 -> Finding the optimal value for beta'

execfile('/home/xb899100/bin/stationParse.py')# <- get wind data
d_method = getDownMethod[experiment]
u_method = getUpMethod[experiment]
differences_NT = []
differences_CT = []
for date in manual_NT:
    if date != '':
        #get the image
        dt_date = dt.strptime(date, '%Y-%m-%d %H:%M:%S')
        image = getOneFrequency(dt_date)
        #get the wind direction myBin
        time = str(dt_date.time()).split(':')
        #get the date in the format that the wind directions are stored in 
        visdate = date + 'UTC'
        #round the minutes to be at the times that there are wind obs
        #convert to a cloud mask
        cloud_mask = image*1.0
        cloud_mask[image > cldmask] = 1.0
        cloud_mask[image <= cldmask] = 0.0
        #round the minutes to be at the times that there are wind obs
        if time[1] == '15':
            visdate_wind = str(dt_date.date()) + ' ' + time[0] + ':20:00'
        else:
            visdate_wind = str(dt_date.date()) + ' ' + time[0] + ':50:00'
        iB = np.where(np.array(myDates) == visdate_wind)[0]
        myBin = myBins[iB]
        difference = get_downwind(method = d_method, cm = cloud_mask, myBin = myBin) - get_upwind(method = u_method, cm = cloud_mask, myBin = myBin)
        differences_NT.append(difference)

for date in manual_CT:
    if date != '':
        #get the image
        dt_date = dt.strptime(date, '%Y-%m-%d %H:%M:%S')
        image = getOneFrequency(dt_date)
        #get the wind direction myBin
        time = str(dt_date.time()).split(':')
        #get the date in the format that the wind directions are stored in 
        visdate = date + 'UTC'
        #round the minutes to be at the times that there are wind obs
        #convert to a cloud mask
        cloud_mask = image*1.0
        cloud_mask[image > cldmask] = 1.0
        cloud_mask[image <= cldmask] = 0.0
        #round the minutes to be at the times that there are wind obs
        if time[1] == '15':
            visdate_wind = str(dt_date.date()) + ' ' + time[0] + ':20:00'
        else:
            visdate_wind = str(dt_date.date()) + ' ' + time[0] + ':50:00'
        iB = np.where(np.array(myDates) == visdate_wind)[0]
        myBin = myBins[iB]
        difference = get_downwind(method = d_method, cm = cloud_mask, myBin = myBin) - get_upwind(method = u_method, cm = cloud_mask, myBin = myBin)
        differences_CT.append(difference)

NT = plt.hist(differences_NT, bins = np.arange(-1.0, 1.01, 0.01), cumulative = -1, normed = True, label = 'NT dist')
CT = plt.hist(differences_CT, bins = np.arange(-1.0, 1.01, 0.01), cumulative = 1, normed = True, alpha = 0.65, label = 'CT dist')
beta = NT[1][np.where(np.abs(NT[0] - CT[0]) == np.nanmin(np.abs(NT[0] - CT[0])))[0][0]]
plt.show()
print 'getAlgorithmThresholds2 -> beta = ' + str(beta)
