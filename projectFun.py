###My Project Functions
#Contains all the functions for my project
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from os import listdir
import datetime


### Send a completion email ###
def send_email(message, subject, attachments, isAttach = True):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.image import MIMEImage
    from email.mime.multipart import MIMEMultipart
    
    # Create the container (outer) email message.
    msg = MIMEMultipart()
    msg['Subject'] = subject
    
    me = 'm.c.johnston@pgr.reading.ac.uk'
    
    msg['From'] = me
    msg['To'] = me
    msg.preamble = '\n'
    
    if isAttach:
        # Assume we know that the image files are all in PNG format
        if type(attachments) != list:
            attachments = [attachments]
        for my_file in attachments:
            # Open the files in binary mode.  Let the MIMEImage class automatically
            # guess the specific image type.
            fp = open(my_file, 'rb')
            img = MIMEImage(fp.read())
            fp.close()
            msg.attach(img)
    
    # Create a text/plain message
    body = MIMEText(message) # convert the body to a MIME compatible string
    msg.attach(body) # attach it to your main message
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(me, me, msg.as_string())
    s.quit()
    
    return None

###Plotting Visible imagery
def getReflectance(date, time, lon, lat, data):
    "takes the calibrated and corrected nominal reflectance values and converts"
    " that to relfectance using the calculated solar zenith angle and sun-earth"
    " distance at each date and time."
    time = int(time)
    year = date/1000 #year of image
    N = date - year*1000 #day of year for that image
    hour = time/10000 #getting the hour and portion of the variable 'time'
    minute = time/100 - hour*100 #getting the minutes portion of the variable 'time'
    hour = hour + minute/60. # getting the decimal hour.
    N = N + hour/24. #getting the decimal day of year with fractions of hour added.
    
    d = 1 + 0.017*np.cos((N - 2)*2*np.pi/365) #perihelion is on the second of Jan 2016
    alpha = (N - 10)*2*np.pi/365
    delta = -23.44*np.cos(alpha)*2*np.pi/360 #approximated solar declination angle

    solarHour = hour - ((720 - 4*lon)/1440)*24 #calculate solar hour by subtracting solar noon
    h = 2*np.pi*solarHour/24 #hour angle
    
    #calculate the cosine of solar zenith angle
    cos_szn = np.sin(lat*2*np.pi/360)*np.sin(delta) + np.cos(lat*2*np.pi/360)*np.cos(delta)*np.cos(h)
    
    #use solar zenith angle to correct data
    data = (data*d*d)/cos_szn
    
    #further limit data to a range 0 to 1.
    data1 = data.copy()
    data1[data1 >= 1.0] = 1.0
    data1[data1 < 0.0] = 0.0
    
    return data1

def getSZA(date, time, lon, lat):
    """
    designed for filtering out the low sun angle images in CTID_vis
    returns the solar zenith angle in degrees
    accepts date strings in the form yyyy-mm-dd
    accepts time strings in the form hh:mm:ss utc
    """
    import datetime
    time0 = time.split(':')
    time = int(time0[0]+time0[1]+time0[2])
    year = date.split('-')[0] #year of image
    date0 = date.split('-')
    N = datetime.datetime(int(date0[0]), int(date0[1]), int(date0[2])).timetuple().tm_yday #day of year for that image
    hour = time/10000 #getting the hour and portion of the variable 'time'
    minute = time/100 - hour*100 #getting the minutes portion of the variable 'time'
    hour = hour + minute/60. # getting the decimal hour.
    N = N + hour/24. #getting the decimal day of year with fractions of hour added.
    
    d = 1 + 0.017*np.cos((N - 2.)*2.*np.pi/365.) #perihelion is on the second of Jan 2016
    alpha = (N - 10)*2*np.pi/365
    delta = -23.44*np.cos(alpha)*2*np.pi/360. #approximated solar declination angle

    solarHour = hour - ((720. - 4.*lon)/1440.)*24. #calculate solar hour by subtracting solar noon
    h = 2.*np.pi*solarHour/24. #hour angle
    
    #calculate the cosine of solar zenith angle
    cos_szn = np.sin(lat*2.*np.pi/360.)*np.sin(delta) + np.cos(lat*2.*np.pi/360.)*np.cos(delta)*np.cos(h)
    SZA = np.arccos(cos_szn)*360./(2.*np.pi)
    return SZA
    
def calibrateVis(data, drift = 1.0):
    "'data' is a netcdf of visible satellite imagery"
    "'drift' is the post-launch calibration for instrument drift"
    #calibrate the data
    #step 1 divide by 32 to convert from 16-bit to 10-bit
    data = data.variables['data'][0]/32
    
    #step 2 apply calibration curve to get nominal reflectance
    Xspace = 29
    k = 0.001160
    data = k*(data - Xspace)
    
    #step 3 post-launch calibration for instrument drift
    data = drift*data
    return data

###Masking functions
def LandMask(elev = -100.):
    DEM = Dataset('/home/xb899100/bin/BermudaDEM.nc')
    DEMlon = DEM.variables['longitude'][:]
    DEMlat = DEM.variables['latitude'][:]
    dem = DEM.variables['elevation'][0,:,:]
    #X,Y = np.meshgrid(DEMlon, DEMlat)
    #x,y = m(X,Y)
    #m.contourf(x, y, dem)
    
    landSeaMask = dem.copy()
    landSeaMask[dem > elev] = 1.
    landSeaMask[dem <= elev] = 0.
    return DEMlon, DEMlat, landSeaMask

def landSea(DEMlon, DEMlat, DEM, maskLon, maskLat, mask):
    # we know that bermuda is roughly at 32.3N and -64.8W
    #we know that a bubble of 0.16 degrees radius excludes all of the island
    #centre point of the island
    centreX = -64.8
    centreY = 32.3
    #all the distances from Bermuda in degrees
    R = np.zeros_like(mask)
    for iX in range(mask.shape[0]):
        for iY in range(mask.shape[1]):
            R[iX, iY] =  np.sqrt((maskLon[iX, iY] - centreX)**2 + (maskLat[iX, iY] - centreY)**2)
    R[R < 0.16] = 1.
    #call the distances that are less than 0.16 degrees away = 1
    
    ixs = [i for i in range(len(DEMlon)) if (DEMlon[i] > centreX - 0.16) & (DEMlon[i] < centreX + 0.16)]
    iys = [i for i in range(len(DEMlat)) if (DEMlat[i] > centreY - 0.16) & (DEMlat[i] < centreY + 0.16)]
    lsm = np.zeros_like(mask)
    for iX in range(mask.shape[0]):
        for iY in range(mask.shape[1]):
            if R[iX, iY] == 1.:
                for ix in ixs:
                    for iy in iys:
                        if DEM[iy,ix] == 1:
                            R2 =  np.sqrt((maskLon[iX, iY] - DEMlon[ix])**2 + (maskLat[iX, iY] - DEMlat[iy])**2)
                            if R2 < 1./60.:
                                lsm[iX, iY] = 1.
    return lsm

def intoBins(vals0, res):
    #vals0 is an array
    vals = vals0.copy()
    #sort into bins
    vals[(0 < vals) & (vals <= res/2.)] = 1 #north
    vals[((360 - res/2.) < vals) & (vals <= 360.)] = 1
    for idir in xrange(1, int(360./res)):
        lower = idir*res - res/2.
        upper = lower + res
        vals[(lower < vals) & (vals <= upper)] = idir + 1
    return vals

def Mask2(radius, lon, lat, land = False, inner = 0.16, res = 10.0):
    ###CIRCULAR MASK
    'creates a mask of wedges radiating from the centre point that are centred'
    'on north, going clockwise, requires a radius in degrees, radius'
    mask = np.zeros_like(lon)
    #define the centre points
    centreX = -64.8
    centreY = 32.3
    
    if land:
        lsm = Dataset('/glusterfs/greybls/users/xb899100/SatData/lsm-sat-coords.nc')
        lsm = lsm.variables['mask'][:,:]
        # land = 1, sea = 0
        #change that so land = 0 and sea = 1.
        lsm -= 1
        lsm = np.abs(lsm)
    #calculate the distance of each point in degrees from the centre
    #calculate the angle away from north of each point
    #sort them into bins 1 to 8 for each direction
    for ix in xrange(lon.shape[0]):
        for iy in xrange(lon.shape[1]):
            R = np.sqrt((lat[ix, iy]-centreY)**2 + (lon[ix, iy]-centreX)**2)
            #find the angle
            if (lat[ix, iy] - centreY) > 0:
                if (lon[ix, iy] - centreX) > 0:
                    angle = np.arcsin(abs(lon[ix, iy] - centreX)/R)*180/np.pi
                elif (lon[ix, iy] - centreX) < 0:
                    angle = 360 - np.arcsin(abs(lon[ix, iy] - centreX)/R)*180/np.pi
            elif (lat[ix, iy] - centreY) < 0:
                if (lon[ix, iy] - centreX) > 0:
                    angle = 180 - np.arccos(abs(lat[ix, iy] - centreY)/R)*180/np.pi
                elif (lon[ix, iy] - centreX) < 0:
                    angle = 180 + np.arccos(abs(lat[ix, iy] - centreY)/R)*180/np.pi
            if angle == 0.0:
                angle = 360.0
            #if within the radii, assign an angle
            if land:
                if(R < radius):
                    mask[ix, iy] = angle
            else:
                if (R < radius) & (R > inner):
                    mask[ix, iy] = angle
                else:
                    mask[ix, iy] = 0
    if land:    
        mask = mask*lsm
        
    #sort into bins
    mask = intoBins(mask, res)
    
    return mask

def Mask2Stagger(radius, lon, lat, land = False, inner = 0.16, res = 10.0):
    ###CIRCULAR MASK
    'creates a mask of wedges radiating from the centre point that are centred'
    'on north+res/2, going clockwise, requires a radius in degrees, radius'
    mask = np.zeros_like(lon)
    #define the centre points
    centreX = -64.8
    centreY = 32.3
    
    if land:
        lsm = Dataset('/glusterfs/greybls/users/xb899100/SatData/lsm-sat-coords.nc')
        lsm = lsm.variables['mask'][:,:]
        # land = 1, sea = 0
        # change that so land = 0 and sea = 1.
        lsm -= 1
        lsm = np.abs(lsm)
    #calculate the distance of each point in degrees from the centre
    #calculate the angle away from north of each point
    #sort them into bins 1 to 8 for each direction
    for ix in xrange(lon.shape[0]):
        for iy in xrange(lon.shape[1]):
            R = np.sqrt((lat[ix, iy]-centreY)**2 + (lon[ix, iy]-centreX)**2)
            #find the angle
            # add res/2 to the angle to stagger the bins
            if (lat[ix, iy] - centreY) > 0:
                if (lon[ix, iy] - centreX) > 0:
                    angle = np.arcsin(abs(lon[ix, iy] - centreX)/R)*180/np.pi + res/2.
                elif (lon[ix, iy] - centreX) < 0:
                    angle = 360 - np.arcsin(abs(lon[ix, iy] - centreX)/R)*180/np.pi + res/2.
            elif (lat[ix, iy] - centreY) < 0:
                if (lon[ix, iy] - centreX) > 0:
                    angle = 180 - np.arccos(abs(lat[ix, iy] - centreY)/R)*180/np.pi + res/2.
                elif (lon[ix, iy] - centreX) < 0:
                    angle = 180 + np.arccos(abs(lat[ix, iy] - centreY)/R)*180/np.pi + res/2.
            if angle == 0.0:
                angle = 360.0
            if angle > 360.0:
                angle -= 360.0
            #if within the radii, assign an angle
            if land:
                if(R < radius):
                    mask[ix, iy] = angle
            else:
                if (R < radius) & (R > inner):
                    mask[ix, iy] = angle
                else:
                    mask[ix, iy] = 0
    if land:    
        mask = mask*lsm
        
    #sort into bins
    mask = intoBins(mask, res)
    
    return mask

def Mask3(radius, lon, lat, res = 10.0, inner_r = 0.0):
    ###BUFFER MASK
    'creates a mask of wedges radiating from the centre point that are centred'
    'on the eight cardinal directions, requires a radius in degrees, radius'
    mask = np.zeros_like(lon)
    mask2 = np.zeros_like(lon)
    #define the centre points
    centreX = -64.8
    centreY = 32.3
    
    lsm = Dataset('/glusterfs/greybls/users/xb899100/SatData/lsm-sat-coords.nc')
    lsm_lon = lsm.variables['x'][:]
    lsm_lat = lsm.variables['y'][:]
    lsm = lsm.variables['mask'][:,:]
    #land points = 1
    #sea points = 0
    
    #find the coordinates of the land points
    x_land, y_land = np.where(lsm == 1)
    #at each land point, mask a circle of radius = radius
    for i in range(len(x_land)):
        #get the coordinates of each point
        x = lsm_lon[x_land[i], y_land[i]]
        y = lsm_lat[x_land[i], y_land[i]]
        #find the distance away from each point
        R = np.sqrt((lat-y)**2 + (lon-x)**2)
        #find all the points that are inside of the circle of radius = radius
        #and set them = 1.0
        mask[R <= radius] = 1.0
        mask2[R <= inner_r] = 1.0
    mask -= mask2
    
    #calculate the angle away from north of each point
    #sort them into bins e.g. 1 to 8 for each direction
    check1 = (lat - centreY)
    check2 = (lon - centreX)
    R = np.sqrt((lat-centreY)**2 + (lon-centreX)**2)
    for ix in xrange(lon.shape[0]):
        for iy in xrange(lon.shape[1]):
            #find the angle
            if check1[ix, iy] > 0:
                if check2[ix, iy] > 0:
                    angle = np.arcsin(abs(lon[ix, iy] - centreX)/R[ix, iy])*180/np.pi
                elif check2[ix, iy] < 0:
                    angle = 360 - np.arcsin(abs(lon[ix, iy] - centreX)/R[ix, iy])*180/np.pi
            elif check1[ix, iy] < 0:
                if check2[ix, iy] > 0:
                    angle = 180 - np.arccos(abs(lat[ix, iy] - centreY)/R[ix, iy])*180/np.pi
                elif check2[ix, iy] < 0:
                    angle = 180 + np.arccos(abs(lat[ix, iy] - centreY)/R[ix, iy])*180/np.pi
            if angle == 0.0:
                angle = 360.0
            
            mask[ix, iy] *= angle
        
    #sort into bins
    mask = intoBins(mask, res)
    mask *= np.abs(lsm-1)
    return mask

def getSectorMeans(mask, masked):
    'gets the mean of each sector in a 2D field masked with Mask2.'
    nSectors = int(np.max(mask))
    
    means = np.zeros(nSectors)
    
    for sector in range(nSectors):
        means[sector] = (np.mean(masked[mask==(sector+1)]))
    
    return means

###Calculating weighted means by pressure
def presWeight(a, b, p, U, DIR, units = 'knots'):
    "calculates pressure weighted average wind speed from a to b."
    "where a is the lower boundary and b is the upper boundary of integration, "
    "p = list of pressure levels with corresponding observations"
    "U = list of wind speeds at those pressure levels"
    "DIR = list of wind directions at those presure levels"
    
    u, v = toComponents(U, DIR)
    
    uMean = 0
    vMean = 0
    DeltaP = 0
    for ip in xrange(len(p)):
        if p[ip] >= a:
            iLow = ip
            #the index of the pressure below or equal to a
        if p[ip] > b:
            iHigh = ip + 1
            #the index of the pressure above or equal to b
        if (p[ip] < a) & (p[ip] > b):
            dp = (p[ip+1]+p[ip])/2 - (p[ip-1]+p[ip])/2
            uMean += u[ip]*dp
            vMean += v[ip]*dp
            DeltaP += dp
            
    uMean += u[iLow-1]*((p[iLow+1] + p[iLow])/2 - a)
    vMean += v[iLow-1]*((p[iLow+1] + p[iLow])/2 - a)
    DeltaP += ((p[iLow+1] + p[iLow])/2 - a)
    
    uMean += u[iHigh]*(b - (p[iHigh] + p[iHigh-1])/2)
    vMean += v[iHigh]*(b - (p[iHigh] + p[iHigh-1])/2)
    DeltaP += (b - (p[iHigh] + p[iHigh-1])/2)
        
    uMean = uMean/(b - a)
    vMean = vMean/(b - a)
    
    UMean, DIRMean = fromComponents(uMean,vMean)
    if units == 'knots':
        UMean = UMean*0.5144 # convert from kn to m/s
    
    return UMean, DIRMean

def toComponents(speed, direction):
    "takes the wind speed and wind direction and breaks the wind speed into"
    "u- and v- components to be returned."
    "speed = wind speed"
    "direction = wind direction"
    u = - speed*np.sin(direction*np.pi/180)
    v = - speed*np.cos(direction*np.pi/180)
    
    return u, v

def fromComponents(u, v, isList = False):
    "Takes the two horizontal wind components and combines them into a wind"
    "speed and direction to be returned."
    "u = the zonal wind component"
    "v = the meridional wind component"
    speed = np.sqrt(np.array(u)**2 + np.array(v)**2)
    if isList:
        direction = np.arcsin(v/speed)*180/np.pi
        for iV in range(len(u)):
            if speed[iV] > 0:
                if (u[iV] >= 0) and (v[iV] > 0):
                    #if u > 0 and v > 0 SW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] > 0) and (v[iV] <= 0):
                    #if u > 0 and v < 0 NW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] <= 0) and (v[iV] > 0):
                    #if u < 0 and v > 0 SE quadrant
                    direction[iV] = 90.0 + direction[iV]
                elif (u[iV] < 0) and (v[iV] <= 0):
                    #if u < 0 and v < 0 NE quadrant
                    direction[iV] = 90.0 + direction[iV]
            else:
                direction[iV] = 0.
    else:
        direction = np.arcsin(v/speed)*180/np.pi
        if speed > 0:
            if (u >= 0) and (v > 0):
                #if u > 0 and v > 0 SW quadrant
                direction = 270.0 - direction
            elif (u > 0) and (v <= 0):
                #if u > 0 and v < 0 NW quadrant
                direction = 270.0 - direction
            elif (u <= 0) and (v > 0):
                #if u < 0 and v > 0 SE quadrant
                direction = 90.0 + direction
            elif (u < 0) and (v <= 0):
                #if u < 0 and v < 0 NE quadrant
                direction = 90.0 + direction
        else:
            direction = 0.
    return speed, direction

#define a function to interpolate
def interpolateZ(PT, Theta, Z, z, day):
    "takes an array of daily data for a month and interpolates and averages it"
    "PT is the empty array for interpolated data"
    "Theta is the uninterpolated data"
    "Z is the array of heights to interpolate to"
    "z is the array of observed heights"
    "day is the day of the month, or column in observations"
    PT[iZ] = PT[iZ] + (Theta[iz, day] - Theta[iz - 1, day])*(Z[iZ] - z[iz - 1, day])/(z[iz, day] - z[iz - 1, day]) + Theta[iz - 1, day]
    return PT
    
def interpolate(aIN, Z, z, error = -9.99990000e+04):
    """
    aIN is the input uninterpolated data
    Z is the grid to interpolate onto
    z is the grid it is currently on
    takes this data and interpolates aIN from z to Z
    if aIN and z are 2D, it assumes that each column is a unique grid to be 
    interpolated to Z
    """
    from scipy import interpolate
    if len(aIN.shape) == 1:
        missing = [i for i,x in enumerate(aIN) if (x == 0) or (x == error)]
        aIN = np.delete(aIN, missing, 0)
        z = np.delete(z, missing, 0)
        aOUT = interpolate.interp1d(z, aIN, fill_value = 'extrapolate')(Z)
    else:
        aOUT = np.zeros((len(Z), aIN.shape[1]))
        for i in range(aIN.shape[1]):
            missing = [j for j,x in enumerate(aIN[:,i]) if (x == 0) or (x == error)]
            aIN_0 = np.delete(aIN[:,i], missing, 0)
            z_0 = np.delete(z[:,i], missing, 0)
            aOUT[:,i] = interpolate.interp1d(z_0, aIN_0, fill_value = 'extrapolate')(Z)
    return aOUT

def interpolateP(aIN, P, p):
    "aIN is the input uninterpolated data"
    "P is the grid to interpolate onto"
    "p is the grid it is currently on"
    "takes this data and interpolates aIN from z to Z"
    "if aIN and p are 2D, it assumes that each column is a unique grid to be "
    "interpolated to P"
    if len(p.shape) > 1:
        ascents = len(p[0,:]) # number of ascents in dataset
        new_levels = len(P)   # number of vertical pressure levels to interpolate onto
        aOUT = np.zeros((new_levels, ascents)) -9.99990000e+04 
                              # set the ouput array to the missing data value
        for ascent in range(ascents):
            aIN_filtered = []
            pIN_filtered = []
            for p_level in range(len(aIN[:,ascent])):
                observation = aIN[p_level, ascent]
                if (observation != 0.0) and (observation != -9.99990000e+04):
                    # if the observation isn't missing
                    aIN_filtered.append(observation)
                    pIN_filtered.append(p[p_level, ascent])
            aIN_filtered = np.array(aIN_filtered) # convert to an array
            pIN_filtered = np.array(pIN_filtered)
            
            ## do the interpolation ##
            
            for interp_p_level in range(new_levels):
                # find the first observed p that is lower than interp_p_level
                if P[interp_p_level] in pIN_filtered:
                    i = np.where(pIN_filtered == P[interp_p_level])[0][0]
                    aOUT[interp_p_level, ascent] = aIN_filtered[i]
                else:
                    i = 0
                    while (P[interp_p_level] < pIN_filtered[i]) and \
                          (i < (len(pIN_filtered) -1)):
                        i += 1
                    
                    M = (aIN_filtered[i] - aIN_filtered[i-1])/(pIN_filtered[i] - pIN_filtered[i-1])
                    X = (P[interp_p_level] - pIN_filtered[i-1])
                    C = aIN_filtered[i-1]
                    aOUT[interp_p_level, ascent] = M*X + C
                    # print aOUT[interp_p_level, ascent]
    else:
        new_levels = len(P)   # number of vertical pressure levels to interpolate onto
        aOUT = np.zeros((new_levels)) -9.99990000e+04 
                              # set the ouput array to the missing data value
        aIN_filtered = []
        pIN_filtered = []
        for p_level in range(len(aIN)):
            observation = aIN[p_level]
            if (observation != 0.0) and (observation != -9.99990000e+04):
                # if the observation isn't missing
                aIN_filtered.append(observation)
                pIN_filtered.append(p[p_level])
        aIN_filtered = np.array(aIN_filtered) # convert to an array
        pIN_filtered = np.array(pIN_filtered)
        
        ## do the interpolation ##
        
        for interp_p_level in range(new_levels):
            # find the first observed p that is lower than interp_p_level
            if P[interp_p_level] in pIN_filtered:
                i = np.where(pIN_filtered == P[interp_p_level])[0][0]
                aOUT[interp_p_level] = aIN_filtered[i]
            else:
                i = 0
                while (P[interp_p_level] < pIN_filtered[i]) and \
                      (i < (len(pIN_filtered) -1)):
                    i += 1
                
                M = (aIN_filtered[i] - aIN_filtered[i-1])/(pIN_filtered[i] - pIN_filtered[i-1])
                X = (P[interp_p_level] - pIN_filtered[i-1])
                C = aIN_filtered[i-1]
                aOUT[interp_p_level] = M*X + C
                # print aOUT[interp_p_level, ascent]
    return aOUT

def annualMean(aIN):
    "aIN is an array for which the mean for each row is calculated and returned"
    "as a 1D array, aOUT"
    nrows = len(aIN[:,0])
    aOUT = np.zeros(nrows)
    for row in range(nrows):
        aOUT[row] = np.mean(aIN[row,:])
    
    return aOUT
    
def temp2theta(TIN, PIN, P0 = 1000.):
    "TIN: input temperatures (K)"
    "PIN: input pressures (hPa)"
    
    # check that a temperature is passed in Kelvin.
    if type(TIN) == np.ndarray:
        if TIN[0] - 273.15 < 0:
            TIN += 273.15
    elif TIN - 273.15 < -100:
        TIN += 273.15
    
    R = 287. #gas constant for dry air
    cp = 1004. #specific heat capacity for dry air
    theta = TIN*(P0/PIN)**(R/cp)
    return theta
    
def myQfromRH(TIN, RHIN, PIN):
    "TIN: input temperatures in Kelvin"
    "RHIN: input relative humidities"
    "PIN: input pressures"
    if np.min(TIN) < 0:
        TIN += 273.15
    RV = 461.5
    RD = 287.05
    E = RD/RV # epsilon
    e0 = 611.14 # Pa
    L = 2.501e06
    PIN = PIN*100. #hPa
    
    es = e0*np.exp((L/RV)*((1./273.16) - (1./TIN)))
    ev = RHIN*es/100.
    q = ev*E/(PIN-ev*(1-E))
    
    return q

def getQ(TIN, RHIN, PIN):
    "TIN: input temperatures (deg C)"
    "RHIN: input relative humidities (%)"
    "PIN: input pressures (hPa)"
    "breaks at temps much lower than -80 C,I assume q = 0 here"
    if (type(TIN) == int) or (type(TIN) in [float, np.float16, np.float32, np.float64]):
        TIN = np.array([TIN])
    if (type(RHIN) == int) or (type(RHIN) in [float, np.float16, np.float32, np.float64]):
        RHIN = np.array([RHIN])
    if (type(PIN) == int) or (type(PIN) in [float, np.float16, np.float32, np.float64]):
        PIN = np.array([PIN])
    #convert to degC
    if np.nanmax(TIN) > 100.:
        TIN = TIN[:] - 273.15
    #convert to pa
    if np.nanmax(PIN) < 10000.:
        p = PIN[:]*100.
    else:
        p = PIN[:]*1.
    #step one - calculate es    
    #set-up the constants
    c0=0.6105851e+03
    c1=0.4440316e+02
    c2=0.1430341e+01
    c3=0.2641412e-01
    c4=0.2995057e-03
    c5=0.2031998e-05
    c6=0.6936113e-08
    c7=0.2564861e-11
    c8=-.3704404e-13
    #calculate the saturation vapor pressure
    es = c0 + TIN*(c1 + TIN*(c2 + TIN*(c3 + TIN*(c4 + TIN*(c5 + TIN*(c6 + TIN*(c7 + TIN*c8)))))))
    
    #step-two calculate the actual vapor pressure
    ev = RHIN*es/100. #in hPa
    
    #step three calculate the specific humidity
    #more constants
    RV = 461.5
    RD = 287.05
    E = RD/RV # epsilon
    q = ev*E/(p-ev*(1-E)) #kg kg-1
    
    i = np.where(TIN < -80.)
    q[i] = 0.
    
    return q

def getQwIce(TIN, RHIN, PIN):
    "TIN: input temperatures (deg C)"
    "QIN: input relative humidities (%)"
    "PIN: input pressures (hPa)"
    "based on teten code and formulation"
    p = PIN*100. #convert to Pa
    temps = TIN + 273.15 # convert to Kelvin
    #constants and thresholds
    e0 = 611.14 # saturation vapor pressure at t0
    t0 = 273.16 # reference temperature in K
    ti = 250.16 # temperature threshold for ice in K
    i1 = 21.875 # dimensionless constant for ice
    i2 = 7.66 # temperature coefficient for ice (K)
    w1 = 17.269 # dimensionless constant for water
    w2 = 35.86 # temperature coefficient for water (k)
    #function to calculate saturation vapor pressure
    def get_es(state, temp):
        if state == 'ice':
            es = e0*(np.exp(i1*(temp - t0)/(temp - i2)))
        elif state == 'water':
            es = e0*(np.exp(w1*(temp - t0)/(temp - w2)))
        
        return es
    myes = np.zeros_like(TIN)
    for ob in range(len(TIN)):
        if temps[ob] <= ti:
            myes[ob] = get_es('ice', temps[ob])
        elif temps[ob] >= t0:
            myes[ob] = get_es('water', temps[ob])
        else:
            myes[ob] = get_es('ice', temps[ob]) + (get_es('water', temps[ob]) - get_es('ice', temps[ob]))*((temps[ob] - ti)/23)**2
    #calculate q
    ev = RHIN*myes/100.
    RV = 461.5
    RD = 287.04
    E = RD/RV # epsilon
    q = ev*E/(p-ev*(1-E)) #kg kg-1
    
    return q

def PTtoTemp(TIN, PIN):
    "converts potential temperature to temperature"
    Rd = 287.05
    cp = 1005.
    Temperature = TIN/((1000./PIN)**(Rd/cp))
    
    return Temperature

def dew2RH(dew, temp, pres):
    "converts a given dew point and temperature at a given pressure to relative humidity"
    #expects inputs in arrauys
    if type(dew) != np.ndarray:
        dew = np.array(dew)
    if type(temp) != np.ndarray:
        temp = np.array(temp)
    if type(pres) != np.ndarray:
        pres = np.array(pres)
    
    #checks if Ts are in Kelvin, if they are - convert to C
    if dew[0]-273.15 > 0:
        dew -= 273.15
    if temp[0]-273.15 > 0:
        temp -= 273.15
    
    q = getQ(dew, 100., pres)
    qs = getQ(temp, 100., pres)
    RH = q/qs
    return RH

def getDew(QIN, PIN1):
    "converts the specific humidity to a vapour pressure"
    "Uses the Bolton equation to calculate the dew point from e"
    PIN = PIN1*100.
    R = 8.314
    Rd = R/0.0289
    Rv = R/0.018
    Lv = 2.501e06
    E = Rd/Rv
    W = QIN/(1. - QIN)
    es0 = 611.2
    
    e = W*PIN/(E + W)
    Td = ((1./273.15) - (Rv/Lv)*np.log(e/es0))**(-1.)
    
    return Td

def getGM(TIN, PIN):
    """
    Calculate the moist adiabatic lapse rate using temperature and pressure
    """
    if (TIN-273.15) < 0:
        T_C = TIN*1.
    else:
        T_C = TIN - 273.15
    QIN  = getQ(T_C, 100., PIN)[0]
    QIN = QIN/(QIN + 1)
    g    = 9.81 #gravity
    Lv   = 2.501e06 #enthalpy of vaporisation
    R    = 287.05 #specific gas constant for dry air
    E    = 0.622
    cpd  = 1005.
    
    A = (Lv*QIN)/(R*TIN)
    B = (E*QIN*Lv**2.)/(R*TIN**2.)
    
    G_m = g*(1. + A)/(cpd + B)
    
    return G_m
    
def getLFC(Tp, Te, QIN, PIN):
    """
    Function to calculate the LFC.
    Tp = the temperature of a parcel, K
    Te = the temperature of the environment, K
    QIN = the specific humidity profile, kg/kg
    PIN = the pressure profile, hPa
    """
    # Locally define some constants
    g    = 9.81 # gravity
    R    = 287.05 # specific gas constant for dry air
    cpd  = 1005. # specific heat capacity of dry air at constant pressure
    G_d  = g/cpd # dry adiabatic lapse rate
    
    # Lets get some unit conversions going
    Q_0  = np.mean(QIN[0]) # surface specific humidity
    T_d0 = getDew(Q_0, PIN[0]) # surface dew point
    
    # Use the University of Wyoming's formula (Bolton, 1980) to find LCL temperature, K
    A    = 1./(T_d0 - 56.)
    B    = np.log(Tp/T_d0)/800.
    LCL  = (1./(A + B)) + 56.
    
    # Begin lifting the parcel to find where it is positively buoyant
    for level in range(len(Te)-1):
        # for each observed vertical level
        p_level_half = 0.5*(PIN[level+1] + PIN[level]) # pressure in between current level and the level above
        rho = (p_level_half*100.)/(R*Tp) # air density in the layer between the current level and the level above
        dz = -(PIN[level+1] - PIN[level])*100./(rho*g) # hydrostatic balance to get the height between the two levels e.g. dp/rho*g, rho = p/RT
        if Tp > LCL:
            #print 'Parcel is below the LCL'
            # if the parcel is warmer than the LCL it is below the LCL and should cool at the dry adiabatic lapse rate
            Tp -= G_d*dz # dry adiabatic lapse rate
            if Tp < LCL:
                # the new parcel is now above the LCL
                # undo the ascent
                Tp += G_d*dz
                # find the difference between the parcel temperature and the LCL temperature
                dT = Tp - LCL
                # find the height change needed for the parcel to get to the LCL temperature
                delta_z = dT/G_d
                # set parcel temperature to the LCL and then lift the parcel the remaining distance but moist adiabatically
                Tp = LCL
                p_level_middle = 0.5*(PIN[level] + (PIN[level] - PIN[level-1])*delta_z/dz + PIN[level-1])
                Tp -= getGM(LCL, p_level_middle)*(dz - delta_z)
        else:
            # above LCL cool at moist lapse rate
            Tp -= getGM(Tp, p_level_half)*dz
            # print [PIN[level], Tp, Te[level], level]
            if Tp >= Te[level+1]:
                # if the parcel is warmer than environment it is free
                LFCT = Tp # set LFC temperature
                break
    if level == (len(Te)-2):
        LFCT = None
    
    return LFCT

def getCAPE(TIN1, QIN, PIN, parcel_type = 0):
    """
    Function to calculate surface based CAPE.
    Tp = User provided parcel (e.g. surface based virtual temperature)
    TIN1 = User provided profile of temperatures in K
    QIN  = User provided profile of specific humidity in kg/kg
    PIN  = User provided profile of pressure in hPa
    """
    from scipy import interpolate
    
    # Convert to the correct units
    if (TIN1[0] - 273.15) < 0:
        Te = TIN1 + 273.15 # converts the environment temperature to Kelvin
    else:
        Te = TIN1*1.
    
    # Locally define some constants
    g = 9.81 # gravity
    cpd = 1005. # specific heat capacity of dry air
    R = 287.05 # gas constant for dry air
    G_d = g/cpd # dry adiabatic lapse rate
    
    ### Parcel Type
    if parcel_type == 0:
        # surface based parcel
        Tp = Te[0]
        Q_0  = QIN[0] # get the surface specific humidity
        T_d0 = getDew(Q_0, PIN[0]) # get the surface dew point
    elif parcel_type == 1:
        # mean over lowest ~ 500 m
        rho = PIN[0]/(R*Te[0])
        dp500 = 500.*rho*g/100.
        p_dummy = PIN[0] - np.arange(0., dp500+0.1)
        Tp = PTtoTemp(np.mean(interpolate.interp1d(PIN, temp2theta(Te, PIN))(p_dummy)), PIN[0])
        Q_0  = np.mean(interpolate.interp1d(PIN, QIN)(p_dummy)) # get the surface specific humidity
        T_d0 = getDew(Q_0, PIN[0]) # get the surface dew point
    else:
        print 'WARNING: Unknown parcel type!!!'
    
    # Calculate the LCL temperature using the University of Wyoming's method (Bolton, 1980)
    A = 1./(T_d0 - 56.)
    B = np.log(Tp/T_d0)/800.
    LCL  = (1./(A + B)) + 56.
    
    ## Find the LFC ##
    LFC  = getLFC(Tp, Te, QIN, PIN) # parcel temperature at level of free convection
    print 'LFC = ' + str(LFC)
    
    ## Initialise the parcel ascent list and the CAPE to accumulate ##
    Tps  = [Tp] # parcel temperature
    CAPE = 0 # convective available potential energy
    CIN  = 0 # Convective Inhibition
    
    ### Start accumulating CAPE ###
    #print LFC
    if LFC != None:
        # an LFC has been found
        for level in range(1,len(Te)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            # convert the thickness of the layer from pressure to height units
            rho = (p_level_half*100.)/(R*Tp)
            dz = - (PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT

            if PIN[level] != 0.0:
                #print round(Tp, 1)
                if Tp >= LCL:
                    #print 'Parcel is below the LCL'
                    # if the parcel is warmer than the LCL it is below the LCL and should cool at the dry adiabatic lapse rate
                    Tp -= G_d*dz # dry adiabatic lapse rate
                    if Tp >= LCL:
                        # if the new parcel temperature is still below the LCL append it and move on
                        Tps.append(Tp)
                    else:
                        # the new parcel is now above the LCL
                        # undo the ascent
                        Tp += G_d*dz
                        # find the difference between the parcel temperature and the LCL temperature
                        dT = Tp - LCL
                        # find the height change needed for the parcel to get to the LCL temperature
                        delta_z = dT/G_d
                        # set parcel temperature to the LCL and then lift the parcel the remaining distance but moist adiabatically
                        Tp = LCL
                        p_level_middle = 0.5*(PIN[level] + (PIN[level] - PIN[level-1])*delta_z/dz + PIN[level-1])
                        Tp -= getGM(LCL, p_level_middle)*(dz - delta_z)
                        # append the correct parcel temperature
                        Tps.append(Tp)
                elif Tp >= LFC:
                    #print 'Parcel is above the LCL and below the LFC'
                    # if the parcel is cooler than the LCL but warmer than the LFC it should cool at the moist adiabatic lapse rate, but not accumulate CAPE
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                    Tps.append(Tp)
                    
                    Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                    Te_level_half = 0.5*(Te[level] + Te[level-1])
                    CIN += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
                else:
                    #print 'Parcel is above the LFC'
                    # otherwise cool at moist adiabatic lapse rate and accumulate CAPE
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                    Tps.append(Tp)
                    
                    # get the mean temperature between levels for the parcel and environment
                    Tp_level_half = 0.5*(Tps[-1] + Tps[-2])
                    Te_level_half = 0.5*(Te[level] + Te[level-1])
                    #print [round(Tp_level_half, 1), round(Te_level_half, 1)]
                    if Tp_level_half > Te_level_half:
                        # if the parcel is warmer than the environment, accumulate CAPE
                        CAPE += g*(Tp_level_half - Te_level_half)*dz/Te_level_half
            else:
                Tp -= G_d*dz
                Tps.append(Tp)
    else:
        # there is no LFC
        for level in range(1,len(Te)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            # convert the thickness of the layer from pressure to height units
            rho = (p_level_half*100.)/(R*Tp)
            dz = - (PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT
            
            if PIN[level] != 0.0:
                if Tp >= LCL:
                    # if below the LCL, cool dry adiabatically
                    Tp -= G_d*dz
                else:
                    # above the LCL and cool moist adiabatcally
                    Tp -= getGM(Tps[-1], p_level_half)*dz
                Tps.append(Tp)
            else:
                Tp -= G_d*dz
                Tps.append(Tp)
    
    # calculate the pressure of the LCL and the LFC
    LCLp = interpolate.interp1d(Tps, PIN, fill_value = 'extrapolate')(LCL)#PIN[0]*(LCL/Te[0])**(cpd/R)
    print LCLp
    LFCp = interpolate.interp1d(Tps, PIN, fill_value = 'extrapolate')(LFC)#PIN[0]*(LFC/Te[0])**(cpd/R)
    print LFCp
    return CAPE, CIN, Tps, LCLp, LFCp

def getCIN(Tvp, TIN1, QIN, PIN):
    "Function to calculate CIN"
    TIN  = TIN1 + 273.15
    g    = 9.81 #gravity
    Lv   = 2.5e06 #enthalpy of vaporisation
    R    = 287 #specific gas constant for dry air
    E    = 0.622
    cpd  = 1004
    G_d  = g/cpd
    Tve  = TIN*(1 + 0.608*QIN) #environmental virtual temperature
    Td   = getDew(QIN, PIN)
    Q_0  = np.mean(Q[0:11])
    T_d0 = getDew(Q_0, PIN[0])
    A = 1./(T_d0 - 56.)
    B = np.log(Tvp/T_d0)/800.
    LCL  = (1./(A + B)) + 56.
    LFC  = getLFC(Tvp, TIN, QIN, PIN) #parcel virtual temperature at level of free convection
    #print LFC
    CIN  = 0 # convective inhibition
    Tvps = [Tvp]
    if LFC != None:
        for level in range(1,len(TIN)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            rho = (p_level_half*100.)/(R*Tvp)
            dz = -(PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT
            if PIN[level] != 0.0:
                if Tvp >= LCL:
                    #if below (warmer than) LCL cool at dry lapse rate
                    Tvp -= G_d*dz
                    Tvps.append(Tvp)
                    #calculate the half points
                    Tvp_level_half = 0.5*(Tvps[-1] + Tvps[-2])
                    Tve_level_half = 0.5*(Tve[level] + Tve[level-1])
                    
                    if Tvp_level_half < Tve_level_half:
                        #if there is negative buoyancy sum CIN
                        CIN += g*(Tvp_level_half - Tve_level_half)*dz/Tve_level_half
                elif Tvp >= LFC:
                    #if above (cooler than) LCL, but below (warmer than) LFC cool
                    #at moist lapse rate
                    Tvp -= getGM(Tvp, p_level_half)*dz
                    Tvps.append(Tvp)
                    #calculate the half points
                    Tvp_level_half = 0.5*(Tvps[-1] + Tvps[-2])
                    Tve_level_half = 0.5*(Tve[level] + Tve[level-1])
                    if Tvp_level_half < Tve_level_half:
                        #if there is negative buoyancy sum CIN
                        CIN += g*(Tvp_level_half - Tve_level_half)*dz/Tve_level_half
                else:
                    #above the LFC and cools at moist lapse rate, stop counting CIN
                    Tvp -= getGM(Tvp, p_level_half)*dz
                    Tvps.append(Tvp)
    else:
        for level in range(1,len(TIN)):
            p_level_half = 0.5*(PIN[level] + PIN[level-1])
            rho = (p_level_half*100.)/(R*Tvp)
            dz = -(PIN[level] - PIN[level-1])*100./(rho*g) # dp/rho*g, rho = p/RT
            if PIN[level] != 0.0:
                if Tvp >= LCL:
                    Tvp -= G_d*dz
                else:
                    Tvp -= getGM(Tvp, p_level_half)*dz
                Tvps.append(Tvp)
            else:
                Tvp -= G_d*dz
                Tvps.append(Tvp)
                        
    return CIN

def getNewCAPE(temp, q_env, pres, height):
    # new method for integrating CAPE
    #check that temperature is in Kelvin
    if (temp[0] - 273.15) < 0:
        temp += 273.15
    
    ## define some constants
    g = 9.81 #gravity
    cpd = 1004. #specific heat capacity of dry air
    Rd = 287. #gas constant for dry air
    
    ## calculate some other neccessities
    #find the index for nearest level to 500m above ground level
    z500 = (500.+height[0])
    i500 = np.where(abs(height - z500) == np.min(abs(height - z500)))[0][0]
    
    q_p0 = 0.0
    for iz in range(i500):
        # assume a linear relationship between levels
        m = (q_env[iz + 1] - q_env[iz])/(height[iz + 1] - height[iz])
        c = q_env[iz]
        if (iz + 1) != i500:
            q_p0 += (0.5*m*height[iz + 1]**2 + c*height[iz + 1]) - (0.5*m*height[iz]**2 + c*height[iz])# average dewpoint of the lowest 500m
        else:
            q_p0 += (0.5*m*z500**2 + c*z500) - (0.5*m*height[iz]**2 + c*height[iz])# average dewpoint of the lowest 500m
    q_p0 /= 500. # divide by 500 to get the average
    td_p0 = getDew(q_p0, pres[0])
    
    
    # calcualate the parcel ascent
    # use the surface (lowest level) based values
    ## 1. calculate the parcel ascent
    theta = temp2theta(temp, pres)
    t_p0 = 0.0
    for iz in range(i500):
        # assume a linear relationship between levels
        m = (theta[iz + 1] - theta[iz])/(height[iz + 1] - height[iz])
        c = theta[iz]
        if (iz + 1) != i500:
            t_p0 += (0.5*m*height[iz + 1]**2 + c*height[iz + 1]) - (0.5*m*height[iz]**2 + c*height[iz])
        else:
            t_p0 += (0.5*m*z500**2 + c*z500) - (0.5*m*height[iz]**2 + c*height[iz]) # average temp of the lowest 500m
    t_p0 /= 500. # divide by 500 to get the average# average temperature of the lowest 500m
    
    # Find the LCL
    A = 1./(td_p0 - 56.)
    B = np.log(t_p0/td_p0)/800.
    LCL  = (1./(A + B)) + 56.
    LCLp = pres[0]*(LCL/temp[0])**(1/0.2854)
    print 'LCLT is ' + str(round(LCL, 1)) + 'K'
    print 'LCLP is ' + str(round(LCLp, 1)) + 'K'
    
    t_par = [temp[0]]
    CAPE = 0
    dZ = 5.
    for level in range(1,len(pres)):
        m = (pres[level] - pres[level - 1])/(height[level] - height[level - 1])
        c = pres[level - 1]
        p_half = (1./(height[level] - height[level - 1]))*(0.5*m*(height[level] - height[level - 1])**2 + c*(height[level] - height[level - 1]))
        if pres[level] > LCLp:
            # use poisson's relationship to estimate the temperature at the new pressure level
            # i.e. cool at dry adiabatic lapse rate below lcl
            t_p = temp[0]/((pres[0]/pres[level])**(Rd/cpd))
            t_par.append(t_p)
        elif (pres[level] <= LCLp) and (pres[level - 1] > LCLp):
            p_half = 0.5*(pres[level] + LCLp)
            rho = p_half*100./(Rd*LCL*(1 + 0.608*getQ(LCL, 100., p_half)[0]))
            dz  = - (pres[level] - LCLp)*100./(g*rho)
            if dz < dZ:
                t_p = t_par[-1] - getGM(t_par[-1], p_half)*dz
            else:
                p_half = pres[level-1]
                t_p = t_par[-1]
                while dz > dZ:
                    p_half = 0.5*(p_half + m*dZ/2. + p_half)
                    t_p = LCL - getGM(LCL, p_half)*dz
                    dz -= dZ
            t_par.append(t_p)
        else:
            # cool at pseudoadiabatic lapse rate above lcl
            rho = p_half*100./(Rd*t_par[-1]*(1 + 0.608*getQ(t_par[-1], 100., p_half)[0]))
            dz = - (pres[level] - pres[level-1])*100./(g*rho)
            if dz < dZ:
                t_p = t_par[-1] - getGM(t_par[-1], p_half)*dz
            else:
                p_half = pres[level-1]
                t_p = t_par[-1]
                while dz > dZ:
                    p_half = 0.5*(p_half + m*dZ/2. + p_half)
                    t_p -= getGM(t_p, p_half)*dZ
                    dz -= dZ
                t_p -= getGM(t_p, p_half)*dz
            t_par.append(t_p)
    ## 2. integrate the positive area from first LFC to last LNB
    # we want to work with respect to virtual potential temperatures
    t_par = np.array(t_par)
    q_par = getQ(t_par*1., 100., pres*1.)
    tv_par = t_par*(1. + 0.608*q_par)
    tv_env = temp*(1. + 0.608*q_env)
    CAPE = 0.0
    for iz in range(len(height)-1):
        # assume a linear relationship between levels
        m1 = (tv_par[iz + 1] - tv_par[iz])/(height[iz + 1] - height[iz])
        c1 = tv_par[iz]
        m2 = (tv_env[iz + 1] - tv_env[iz])/(height[iz + 1] - height[iz])
        c2 = tv_env[iz]
        if pres[iz] <= LCLp:
            if (tv_par[iz] >= tv_env[iz]) and (tv_par[iz + 1] >= tv_env[iz + 1]):
                intB = 0.5*(m1 - m2)*(height[iz + 1] - height[iz])**2 + (c1 - c2)*(height[iz + 1] - height[iz])
                tv_env_half = (1./(height[iz + 1] - height[iz]))*(0.5*m2*(height[iz + 1] - height[iz])**2 + c2*(height[iz + 1] - height[iz]))
            elif (tv_par[iz] >= tv_env[iz]) and (tv_par[iz + 1] < tv_env[iz + 1]):
                zNeg = (c2 - c1)/(m1 - m2) + height[iz]
                intB = 0.5*(m1 - m2)*(zNeg - height[iz])**2       + (c1 - c2)*(zNeg - height[iz])
                tv_env_half = (1./(zNeg - height[iz]))*(0.5*m2*(zNeg - height[iz])**2 + c2*(zNeg - height[iz]))
            elif (tv_par[iz] < tv_env[iz]) and (tv_par[iz + 1] >= tv_env[iz + 1]):
                zNeg = (c2 - c1)/(m1 - m2) + height[iz]
                intB = 0.5*(m1 - m2)*(height[iz + 1] - zNeg)**2 + (c1 - c2)*(height[iz + 1] - zNeg)
                tv_env_half = (1./(height[iz + 1] - zNeg))*(0.5*m2*(height[iz + 1] - zNeg)**2 + c2*(height[iz + 1] - zNeg))
            else:
                intB = 0
                tv_env_half = 1.
            CAPE += g*intB/tv_env_half
    print 'CAPEV is ' + str(round(CAPE, 2)) + 'J/kg'
    return CAPE, t_par

# http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature
    
def getCCL(temp, pres, dew = [None], q = [None]):
    """
    Function to find the convective condensation level (CCL)
    CCL is the level at which  the surface dew point equals the environmental 
    temperature. The convective temperature is then the temperature a parcel
    taken from that level would have if brought dry adiabatically to the 
    surface.
    """
    Rd = 287.
    cp = 1004.
    # find the level at which the surface specific humidity equals the 
    # environmental saturation specific humidity
    
    #check that temp is in Kelvin
    if temp[0]-273.15 < 0:
        temp += 273.15
    q_env_sat = getQ(temp, 100., pres)
    if dew[0] != None:
        q = getQ(dew, 100., pres)
    
    i = np.where(abs(q[0] - q_env_sat) == np.min(abs(q[0] - q_env_sat)))[0][0]
    
    # if the temperature at this level is greater than the dew point, then 
    # it is below the CCL and we need to interpolate above, else the opposite
    if q_env_sat[i] > q[0]:
        CCL_p = (q[0] - q_env_sat[i])*(pres[i+1] - pres[i])/(q_env_sat[i+1] - q_env_sat[i]) + pres[i]
    elif q_env_sat[i] < q[0]:
        CCL_p = (q[0] - q_env_sat[i-1])*(pres[i] - pres[i-1])/(q_env_sat[i] - q_env_sat[i-1]) + pres[i-1]
    else:
        CCL_p = pres[i]
    
    CCL_T = getDew(q[0], CCL_p)*(pres[0]/CCL_p)**(Rd/cp)
    return CCL_p, CCL_T
    
def get_LCL_p(temp, dew, pres):
    #these are just the surface temperatures in K
    #just the surface pressurepressure in hPa
    cp = 1004.
    Rd = 287.
    A = 1./(dew - 56.)
    B = np.log(temp/dew)/800.
    LCL_temp  = (1./(A + B)) + 56.
    LCL_pres = pres*(LCL_temp/temp)**(cp/Rd)
    return LCL_pres

def get_ML_p(theta_v, pres):
    #temperature in K
    #pres in hPa
    #uses parcel method
    k = 1
    while (theta_v[k] < theta_v[0]) and (k < len(theta_v)):
        k += 1
    if theta_v[k] > theta_v[0]:
        ML_pres = (pres[k] - pres[k-1])*(theta_v[0] - theta_v[k-1])/(theta_v[k] - theta_v[k-1]) + pres[k-1]
    elif theta_v[k] == theta_v[0]:
        ML_pres = pres[k]
    else:
        k = 1
        while (theta_v[k] < theta_v[1]) and (k < len(theta_v)):
            k += 1
        if theta_v[k] > theta_v[1]:
            ML_pres = (pres[k] - pres[k-1])*(theta_v[1] - theta_v[k-1])/(theta_v[k] - theta_v[k-1]) + pres[k-1]
        elif theta_v[k] == theta_v[1]:
            ML_pres = pres[k]
    
    return ML_pres

def great_circle_distance(lat1, lon1, lat2, lon2, units = 'km'):
    """ Accepts the latitude and longitude of two locations on Earth and 
    calculates the distance between those two points using a Great Circle
    Distance equation. Latitude is provided in degrees North, and Longitude in
    degrees East. """
    import numpy as np
    
    # Convert from degrees to radians
    lat1_rad = lat1*np.pi/180.
    lon1_rad = lon1*np.pi/180.
    lat2_rad = lat2*np.pi/180.
    lon2_rad = lon2*np.pi/180.
    
    # Calculate the chord length
    dX = np.cos(lat2_rad)*np.cos(lon2_rad) - np.cos(lat1_rad)*np.cos(lon1_rad)
    dY = np.cos(lat2_rad)*np.sin(lon2_rad) - np.cos(lat1_rad)*np.sin(lon1_rad)
    dZ = np.sin(lat2_rad) - np.sin(lat1_rad)
    C = np.sqrt(dX**2 + dY**2 + dZ**2)
    
    # Calculate the central angle
    dS = 2.*np.arcsin(C/2.)
    
    # Calculate the Great Circle Distance
    r = 6371. # Approximate mean radius of Earth in km
    d = r*dS
    
    # Convert to different units
    if units in ['km', 'kilometers', 'kilometres', 'k', 'KM', 'Kilometres', 'Kilometers', 'K']:
        return d
    elif units in ['m', 'metres', 'meters', 'met', 'M', 'Metres', 'Meters', 'Met']:
        d *= 1000.
        return d
    elif units in ['sm', 'mi', 'mile', 'miles', 'SM', 'MI', 'Mile', 'Miles']:
        d /= 1.609
        return d
    elif units in ['nm', 'nautical miles', 'nautical mile', 'NM', 'Nautical Miles', 'Nautical Mile']:
        d /= 1.85035
        return d

print 'Reading, UK is about ' + str(round(great_circle_distance(32.3, -64.8, 51.5, -0.95), 0)) + 'km from Bermuda'

