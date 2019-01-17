#get the raw sounding data
import numpy as np
import matplotlib.pyplot as plt
import csv
execfile('/home/xb899100/bin/projectFun.py')

def get_raw_soundings():
    "reads the raw sounding data and removes missing values"
    
    path = '/home/users/xb899100/OLD/Sounding/IGRA/'
    
    #MLZp is the pressure of the top of the mixed layer found using the parcel method
    #LCLp is the pressure of the lifting condensation level
    
    igra_data = []
    with open(path + 'derived/derived.txt', 'r') as igra_file:
        for line in igra_file:
            igra_data.append(line)
    
    
    ndays = 2383 # number of ascents in 2012-2016 period
    soundings = {'height'  : np.zeros((200, ndays)), 
                'pres'     : np.zeros((200, ndays)), 
                'temp'     : np.zeros((200, ndays)),
                'theta'    : np.zeros((200, ndays)),
                'RH'       : np.zeros((200, ndays)), 
                'q'        : np.zeros((200, ndays)),
                'temp_d'   : np.zeros((200, ndays)),
                'uwind'    : np.zeros((200, ndays)), 
                'vwind'    : np.zeros((200, ndays)),
                'wind_spd' : np.zeros((200, ndays)),
                'wind_dir' : np.zeros((200, ndays)),
                'date'     : np.zeros((1, ndays)),
                'INVTD'    : np.zeros((1, ndays)),
                'MLZp'     : np.zeros((1, ndays)),
                'LCLp'     : np.zeros((1, ndays)),
                'CAPE'     : np.zeros((1, ndays)),
                'CIN'      : np.zeros((1, ndays))}
    it = -1
    for line in igra_data:
        myLine = line.split()
        if line[0] == '#' and int(myLine[1]) >= 2012 and int(myLine[1]) <= 2016:
            it += 1
            
            level = igra_data.index(line) + 1
            myInc = 0
            nlevels = int(myLine[6])
            
            soundings['date'][0,it]  = int(myLine[1] + myLine[2] + myLine[3] + myLine[4])
            soundings['INVTD'][0,it] = int(line[56:61])/10.
            soundings['MLZp'][0,it]  = int(line[62:67])/100.#convert to hPa
            soundings['LCLp'][0,it]  = int(line[86:91])/100. #convert to hPa
            soundings['CAPE'][0,it]  = int(line[146:151])
            soundings['CIN'][0,it]   = int(line[152:157])
            for inc in range(nlevels):
                data_line = igra_data[level+inc]
                z = int(data_line[17:23]) #m, the calculated geopotential height
                # QC pressure
                p = int(data_line[1:7])/100. #convert to hPa
                # QC temperatures
                T = int(data_line[25:31])/10. #convert to K
                # QC RH
                RH = int(data_line[97:103])/10. #convert to %
                if RH < 0:
                    RH = -9.99990000e+04
                #QC wind
                u = int(data_line[113:119])/10. #convert to m/s
                v = int(data_line[129:135])/10. #convert to m/s
                if (u < -999.) or (v < -999.):
                    u = -9.99990000e+04
                    v = -9.99990000e+04
                
                #only include complete sounding levels
                if z > 0.:
                    if p > 0.:
                        soundings['height'][myInc, it]   = z
                        soundings['pres'][myInc, it]     = p
                        soundings['temp'][myInc, it]     = T
                        soundings['theta'][myInc, it]    = temp2theta(T,p)
                        soundings['RH'][myInc, it]       = RH
                        if RH != -9.99990000e+04:
                            soundings['q'][myInc, it]    = getQ(T, RH, p)[0]
                            soundings['temp_d'][myInc, it] = getDew(soundings['q'][myInc, it], p)
                        else:
                            soundings['q'][myInc, it]    = -9.99990000e+04
                            soundings['temp_d'][myInc, it] = -9.99990000+04
                        
                        wind = fromComponents(u,v)
                        if (u == -9.99990000e+04) or (v == -9.99990000e+04):
                            wind = [-9.99990000e+04, -9.99990000e+04]
                        soundings['uwind'][myInc, it]    = u
                        soundings['vwind'][myInc, it]    = v
                        soundings['wind_spd'][myInc, it] = wind[0]
                        soundings['wind_dir'][myInc, it] = wind[1]
                        myInc += 1
    return soundings

def get_interpolated_soundings(P, soundings = None):
    "Uses interpolateP() function from project fun to interpolate onto a common grid"
    
    if type(soundings) == 'dict':
        ## Read in the raw soundings ##
        soundings = get_raw_soundings()
    
    ## Create a dictionary for the interpolated soundings"
    interpolated = {}
    
    ## Do the interpolation ##
    for key in soundings.keys():
        if key != 'pres':
            if soundings[key].shape[0] == 200:
                print 'Interpolating ' + key
                interpolated[key] = interpolateP(soundings[key], P, soundings['pres'])

    ## Add the dates
    interpolated['date'] = soundings['date']
    
    return interpolated