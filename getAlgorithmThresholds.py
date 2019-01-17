#Script to generate values for alpha ane beta for use in the algorithm
execfile('/home/xb899100/bin/projectFun.py') #main project functions
execfile('/glusterfs/greybls/users/xb899100/SatData/readReflectances.py')
#function to read satellite imagery data

#Alpha - the mean cloud fraction to check if it is obscured
def get_alpha(method, cm):
    "methods for getting the value of alpha"
    "method 1 = use the whole domain"
    "method 2 = use the mask area around Bermuda"
    if method == 1:
        a = np.nanmean(cm)
    elif method == 2:
        a = np.nanmean(cm*mask/mask)
    
    return a        
#Beta - the downwind minus upwind difference to check if there is a cloud trail
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

#this requires thesectors that are used in the algorithm
lsmQ = True
innerR = 0.16
outerR = 0.25
myres = 10.0
numberRes = int(360./myres)
window = 90/2.

#mask = Mask2(radius = outerR, lon = lon, lat = lat, land = lsmQ, inner = innerR, res = myres)
mask = Mask3(radius = outerR, lon = lon, lat = lat, res = myres)

cloud_fractions = [] #initialise an empty list to store the cloud fractions
differences = [] #initialise an empty list to store the downwind-upwind differences
count = 0
for date in visdates:
    #check that the scene is well lit
    #get the solar zenith angle
    sza_date = str(date.date())
    sza_time = str(date).split()[1]
    SZA = getSZA(sza_date, sza_time, lon, lat)
    if np.max(SZA) <= 75.0:
        #read that date/time's visible image
        visible_data = getOneFrequency(date)
        
        #create a cloud mask
        cloud_mask = visible_data*1.0
        cloud_mask[visible_data >= 0.15] = 1.0
        cloud_mask[visible_data < 0.15] = 0.0
        
        #get that date/time's cloud fraction
        cloud_fraction = get_alpha(method = getAlphaMethod[experiment], cm = cloud_mask)
        #append it to the list
        cloud_fractions.append(cloud_fraction)
        #get all possible upwind/downwind pairs
        for myBin in range(1, numberRes+1):
            #get the upwind
            upwind = get_upwind(method = getUpMethod[experiment], cm = cloud_mask, myBin = myBin)
            #get downwind
            downwind = get_downwind(method = getDownMethod[experiment], cm = cloud_mask, myBin = myBin)
            
            #get the difference
            difference = downwind - upwind
            #append to the list
            differences.append(difference)
        count += 1
#once all values for all images are calculated:
alpha = np.nanmean(cloud_fractions)
beta = np.std(differences)
###WRITE DATA###
import csv
#ave to base_path
with open(base_path + 'cloud_fractions.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter = ',')
    spamwriter.writerow(cloud_fractions)
with open(base_path + 'differences.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter = ',')
    spamwriter.writerow(differences)
####Code for getting the 90 degree downwind window investigated
# #get the downwind direction
# downwinds = []
# adjustment = int(window/numberRes)
# for adj in range(-adjustment, 1+adjustment):
    # down = up - numberRes/2 + adj
    # #prevent direction folding
    # if down < 1:
        # down += numberRes
    # elif down > numberRes:
        # down -= numberRes
    # #get the cloud fraction in the downwind direction
    # downwind = np.nanmean(cloud_mask[mask == down])
    # downwinds.append(downwind)

# downwind = np.nanmax(downwinds)