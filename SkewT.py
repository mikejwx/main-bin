import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
execfile('/home/xb899100/bin/projectFun.py')
#contains methods to plot a skew T when passed temperature and specific humidity values
#varying with height and pressure
#data at altitudes lower than the 100hPa level is plotted

### Constants
g = 9.81 #gravity
cpd = 1004. #specific heat capacity of dry air
cl = 4200. #specific heat capacity of liquid water
R = 8.314 #gas constant
Mrd = 0.0289 #molar mass of dry air
Mrv = 0.018 #molar mass of water vapour
EPS = Mrv/Mrd #ratio of molar masses of water vapour and dry air
Rd = R/Mrd #gas constant for dry air
Rv = R/Mrv #gas constant for water vapour
Lv = 2.5e06 #latent heat of vapourisation
Lm = 3.3e05 #latent heat of melting
Ls = Lv + Lm #latent heat of sublimation
gamma_d = -g/cpd #dry adiabatic lapse rate
p0 = 1000. #standard pressure level
adiabat_p = np.arange(1050., 99., -5.) #pressure levels to calculate adiabats
T0 = np.arange(-40, 300., 20.)
#p0 level temperatures to calculate dry adiabats
T1 = np.arange(-40., 90., 3)

def dry_adiabats():
	#uses T0 to draw dry adiabats along
	#uses pressure levels every 5hPa and poisson equation for potential temperatue
	adiabats = np.zeros((len(adiabat_p), len(T0)))
	for T in range(len(T0)):
		adiabats[:,T] = (T0[T]+273.15)/(p0/adiabat_p)**(Rd/cpd) - 273.15
			
	return adiabats

def get_saturation_mixing_ratio(temp):
	#calculates the mixing ratio
	#requires temperature in deg C
	es = 6.112*np.exp(17.67*temp/(temp+243.5))
	rs = EPS*es/(adiabat_p - es)
	
	return rs

def moist_adiabats():
    #uses T1 to draw moist adiabats along
    #uses pressure levels every 5hPa to calculate an estimate of the equivalent potential temperature
    #using the moist adiabatic lapse rate
    moist_adiabats = np.zeros((len(adiabat_p), len(T1)))
    moist_adiabats[0,:] = T1
    for T in range(1,len(adiabat_p)):
        # loop over the pressure levels
        for t in range(len(T1)):
            rho = (0.5*(adiabat_p[T-1]+adiabat_p[T])*100.)/(Rd*(moist_adiabats[T-1,t] + 273.15)*(1 + 0.608*getQ(moist_adiabats[T-1,t], 100., 0.5*(adiabat_p[T-1]+adiabat_p[T]))[0]))
            dz = (500./(rho*g)) # dp/rho*g, rho = p/RT
            moist_adiabats[T,t] = moist_adiabats[T-1, t] - \
                getGM(moist_adiabats[T-1,t]+273.15, 0.5*(adiabat_p[T-1] + adiabat_p[T]))*dz
    return moist_adiabats

wmo_p = [1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100]

def skew(temperature, pressure):
    #skews the temperatures by Gamma_d/2
    skewedT = np.zeros_like(temperature) + temperature
    for p_level in range(len(pressure)):
        dp = np.log10(pressure[p_level]) - np.log10(1000.)
        skewedT[p_level] -= 5*dp*9.81 # skew 45 degrees is -5 degrees per dp, and the multiply by gravity to make it proportional in log space
    return skewedT


def paper():
    theta = dry_adiabats()
    thetae = moist_adiabats()
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
        
    therms = np.arange(-100., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    fig = plt.gcf()
    fig.set_size_inches(6, 8)
    plt.semilogy(theta, adiabat_p, color = 'g', linewidth = 0.25)
    plt.semilogy(thetae, adiabat_p, color = 'b', linewidth = 0.25)
    plt.semilogy(isotherms, adiabat_p, color = 'gray', linewidth = 0.25, linestyle = '--')
    plt.gca().invert_yaxis()
    plt.ylim([1005, 100])
    plt.gca().yaxis.grid()
    plt.yticks(wmo_p, wmo_p)
    plt.xlim([-40, 40])


def windBarbs(u_in, v_in, p_in, gs):
    #I want to only use winds on the wmo levels
    u_wmo = []
    v_wmo = []
    for i in range(len(p_in)):
        if p_in[i] in wmo_p:
            u_wmo.append(u_in[i])
            v_wmo.append(v_in[i])
    u_wmo = np.array(u_wmo)
    v_wmo = np.array(v_wmo)
    #adds wind barbs to the SkewT
    ax2 = plt.subplot(gs[1])
    ax2.barbs(np.zeros_like(wmo_p), wmo_p, u_wmo, v_wmo, barbcolor = 'k', clip_on=False, zorder=100)
    plt.xlim([-1, 1])
    plt.ylim([1005, 100])
    plt.yticks([])
    plt.xticks([])
    plt.axis('off')

def plotSkewT(temp, t_dew, p, u = np.array([-999]), v = np.array([-999]), 
              parcel_temp = np.array([-999]), CAPE = -999, date = -999, 
              my_title = '', temp_col = ['r'], dew_col = ['g']):
    # plots a skewT on the paper I've defined above.
    # includes a parcel ascent if available
    # includes CAPE in the title if provided
    # temp and t_dew in degC
    # p in hPa
    ### get the correct units for temp and t_dew ###
    # should be in deg C
    if np.min(temp) > 0:
        temp -= 273.15
    if np.min(t_dew) > 0:
        t_dew -= 273.15
    
    gs = gridspec.GridSpec(1, 2, width_ratios = [4, 1])
    gs.update(wspace=0.0, hspace=0.0)
    #paper()
    
    ## plot the skew T paper ##
    theta = dry_adiabats()
    thetae = moist_adiabats()
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
    
    therms = np.arange(-100., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    ax = plt.subplot(gs[0])
    ax.semilogy(theta, adiabat_p, color = 'g', linewidth = 0.25)
    ax.semilogy(thetae, adiabat_p, color = 'b', linewidth = 0.25)
    ax.semilogy(isotherms, adiabat_p, color = 'gray', linewidth = 0.25, linestyle = '--')
    
    ## plot the sounding ##
    if type(temp) != list:
        print 'WARNING! converting to a list of temperature profiles!'
        temp = [temp]
    if type(t_dew) != list:
        print 'WARNING! Converting to a list of dewpoint profiles.'
        t_dew = [t_dew]
    if type(u) != list:
        u = [u]
    if type(v) != list:
        v = [v]
    if len(temp) != len(t_dew):
        print 'You have an unequal number of temperature and dewpoint ascents'
    if len(temp_col) != len(temp):
        print 'WARNING! The number of temperature colors provided does not match the number of temperature profiles. Taking the same dewpoint color for each profile.'
        while len(temp_col) != len(temp):
            temp_col.append(temp_col[-1])
    if len(dew_col) != len(t_dew):
        print 'WARNING! The number of dewpoint colors provided does not match the number of dewpoint profiles. Taking the same dewpoint color for each profile.'
        while len(dew_col) != len(t_dew):
            dew_col.append(dew_col[-1])
    #print('Plotting Temperature')
    
    for i in range(len(temp)):
        ax.plot(skew(temp[i], p), p, color = temp_col[i], lw = 2)
    #print('Plotting Dew point')
        ax.plot(skew(t_dew[i], p), p, color = dew_col[i], ls = '--', lw = 2)
    if parcel_temp[0] != -999:
        #print('Plotting Parcel Ascent')
        if np.min(parcel_temp) > 0:
            parcel_temp -= 273.15
        ax.plot(skew(parcel_temp, p), p, color = 'gray', lw = 2)
    if CAPE != -999:
        #print('Adding CAPE in the title')
        my_title = 'CAPE using $T_v$ = ' + str(round(CAPE, 0)) + '$Jkg^{-1}$'
    
    plt.gca().invert_yaxis()
    plt.ylim([1025, 100])
    plt.gca().yaxis.grid()
    plt.yticks(wmo_p, wmo_p)
    plt.xlim([-40, 40])
    plt.xlabel('Temperature ($^o$C)')
    plt.ylabel('Pressure (hPa)')
    
    if date != -999:
        my_title = date + '\n' + my_title
    plt.title(my_title)    
    ## plot the wind barbs ##
    
    for i in range(len(u)):
        if u[i][0] != -999:
            # windBarbs(u, v, p, gs)
            # only use wind obs on wmo pressure levels so we can see all the wind barbs
            u_wmo = []
            v_wmo = []
            # get the indices corresponding tothe wmo pressure levels to plot the winds
            wmo_indices = []
            for wmo_pres in wmo_p:
                wmo_indices.append(np.where(np.abs(p-wmo_pres) == np.min(np.abs(p-wmo_pres)))[0][0])
            u_wmo = np.array([u[i][index] for index in wmo_indices])
            v_wmo = np.array([v[i][index] for index in wmo_indices])
            # add wind barbs to the SkewT
            ax2 = plt.subplot(gs[1])
            ax2.barbs(np.zeros_like(wmo_p), wmo_p, u_wmo, v_wmo, barbcolor = 'k', \
                clip_on = False, zorder=100)
            plt.xlim([-1, 1])
            plt.ylim([1025, 100])
            print 'Plotted wind barbs'
            plt.yticks([])
            plt.xticks([])
            plt.axis('off')


    
###real data to plot
#path = '/home/xb899100/bin/'
#import scipy.io as io
#data = io.readsav('/home/xb899100/bin/profilearmsonde_save_20010401-20060816.dat')
#sounding = 1130
#temp = data['tdry'][sounding,3:-1]
#Q = data['qsat'][sounding,3:-1]*data['rh'][sounding,3:-1]/100.
#q_kgkg = Q/1000.
#temp_v = temp #*(1 + 0.608*q_kgkg)
#Z = data['alt'][sounding,3:-1]
#p = data['pres'][3:-1]
#t_dew = getDew(q_kgkg, p) - 273.15
#Tvp = np.nanmean(theta(temp_v[0:11], p[0:11]))
#myCAPE = getCAPE(Tvp, temp_v-273.15, q_kgkg, p, Z)
# test_data = []
# with open(path + 'test_sounding.csv', 'r') as test_file:
    # for line in test_file:
        # test_data.append(line)

# p = []
# temp = []
# t_dew = []
# Z = []
# for line in test_data:
	# myline = line.split(',')
	# p.append(float(myline[0]))
	# temp.append(float(myline[2]))
	# t_dew.append(float(myline[3]))
	# Z.append(float(myline[1]))

# execfile('/glusterfs/greybls/users/xb899100/SatData/readCloudFrequencies.py')
# ###Key inputs
# lsmQ = True
# innerR = 0.16
# outerR = 0.25
# myres = 22.5
# lowerP = 950
# upperP = 850
# execfile('/glusterfs/greybls/users/xb899100/SatData/soundingParse.py')
# i = 0
# i += 10
# #p = np.array(p)
# p = interpolated['p'][:,i]
# #temp = np.array(temp)
# temp = interpolated['Temp'][:,i] -273.15
# #t_dew = np.array(t_dew)
# RH = interpolated['RH'][:,i]
# #Z = np.array(Z)


# #get parcel
# #q_kgkg = getQ(t_dew, 100., p)
# q_kgkg = getQ(temp, RH, p)
# t_dew = getDew(q_kgkg, p) - 273.15
# temp_v = (temp+273.15)*(1 + 0.608*q_kgkg)


# CAPE = getCAPE(temp_v[0], temp_v -273.15, q_kgkg, p, Z)

# #CAPE = getCAPE(temp_v[0], temp, q_kgkg, p)

# parcel_temp = np.array(CAPE[1])-273.15
#paper()
#plt.plot(skew(np.array(temp)-273.15, p), p, color = 'r')
#plt.plot(skew(t_dew, p), p, color = 'g')
#plt.plot(skew(np.array(myCAPE[1])-273.15, p), p, color = 'gray')
#plt.title('CAPE using $T$ = ' + str(round(myCAPE[0], 1)) + '$Jkg^{-1}$')
#plt.xlabel('Temperature ($^o$C)')
#plt.ylabel('Pressure (hPa)')
#plt.show()
# plt.savefig('/home/xb899100/bin/test_skewt.png', bbox_tight = True)

