import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import interpolate

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

def getGM(TIN, PIN):
    if (TIN-273.15) < 0:
        T_C = TIN*1.
    else:
        T_C = TIN - 273.15
    QIN  = getQ(T_C, 100., PIN)[0]
    g    = 9.81 #gravity
    Lv   = 2.5e06 #enthalpy of vaporisation
    R    = 287. #specific gas constant for dry air
    E    = 0.622
    cpd  = 1004.
    
    A = (Lv*QIN)/(R*TIN)
    B = (E*QIN*Lv**2)/(R*TIN**2)
    
    G_m = g*(1 + A)/(cpd + B)
    
    return G_m

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

#contains methods to plot a skew T when passed temperature and specific humidity values
#varying with height and pressure
#data at altitudes lower than the 100hPa level is plotted

### Constants
g = 9.81 #gravity
cpd = 1005. #specific heat capacity of dry air
cl = 4200. #specific heat capacity of liquid water
R = 8.314 #gas constant
Mrd = 0.0289 #molar mass of dry air
Mrv = 0.018 #molar mass of water vapour
EPS = Mrv/Mrd #ratio of molar masses of water vapour and dry air
Rd = R/Mrd #gas constant for dry air
Rv = R/Mrv #gas constant for water vapour
Lv = 2.501e06 #latent heat of vapourisation
Lm = 3.3e05 #latent heat of melting
Ls = Lv + Lm #latent heat of sublimation
gamma_d = -g/cpd #dry adiabatic lapse rate
p0 = 1000. #standard pressure level
adiabat_p = np.arange(1050., 99., -5.) #pressure levels to calculate adiabats
w_p = np.arange(1050., 500., -5.) # pressure level on which to plot lines of constant w
target_w = np.array([0.1, 1, 2, 4, 7, 10, 15, 20, 30])*1e-3 # constant mixing ratios to plot lines
T0 = np.arange(-40, 300., 20.)
#p0 level temperatures to calculate dry adiabats
T1 = np.arange(-40., 90.1, 5)

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
    i_1000 = np.where(adiabat_p == 1000)[0][0]
    moist_adiabats[i_1000,:] = T1
    for T in range(i_1000+1,len(adiabat_p)):
        # loop over the pressure levels
        for t in range(len(T1)):
            rho = (0.5*(adiabat_p[T-1]+adiabat_p[T])*100.)/(Rd*(moist_adiabats[T-1,t] + 273.15)*(1 + 0.608*getQ(moist_adiabats[T-1,t], 100., 0.5*(adiabat_p[T-1]+adiabat_p[T]))[0]))
            dz = (500./(rho*g)) # dp/rho*g, rho = p/RT
            moist_adiabats[T,t] = moist_adiabats[T-1, t] - \
                getGM(moist_adiabats[T-1,t]+273.15, 0.5*(adiabat_p[T-1] + adiabat_p[T]))*dz
    for T in range(i_1000, 0, -1):
        # loop over the pressure levels
        for t in range(len(T1)):
            rho = (0.5*(adiabat_p[T-1]+adiabat_p[T])*100.)/(Rd*(moist_adiabats[T,t] + 273.15)*(1 + 0.608*getQ(moist_adiabats[T,t], 100., 0.5*(adiabat_p[T-1]+adiabat_p[T]))[0]))
            dz = (500./(rho*g)) # dp/rho*g, rho = p/RT
            moist_adiabats[T-1,t] = moist_adiabats[T, t] + \
                getGM(moist_adiabats[T,t]+273.15, 0.5*(adiabat_p[T-1] + adiabat_p[T]))*dz

    return moist_adiabats

def const_ws():
    """
    Creates curves for lines with constant specific humidity
    """
    const_ws = np.zeros((len(w_p), len(target_w)))
    for iw in range(len(target_w)):
        const_ws[:,iw] = getDew(target_w[iw]/(1. + target_w[iw]), w_p)-273.15
    
    return const_ws

wmo_p = [1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100]

def skew(temperature, pressure):
    #skews the temperatures by Gamma_d/2
    c = 0.01219513#np.tan(30*np.pi/180.)
    skewedT = np.zeros_like(temperature) + temperature
    for p_i in range(len(pressure)):
        dp = np.log10(pressure[p_i]/1000.)
        skewedT[p_i] -= dp/c
    return skewedT

def paper():
    theta = dry_adiabats()
    thetae = moist_adiabats()
    my_w = const_ws()
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
    for sfc in range(len(target_w)):
        my_w[:,sfc] = skew(my_w[:,sfc], w_p)
    
    therms = np.arange(-120., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    fig = plt.gcf()
    fig.set_size_inches(6, 8)
    plt.semilogy(theta, adiabat_p, color = 'g', lw = 0.25)
    plt.semilogy(thetae, adiabat_p, color = 'b', lw = 0.25)
    plt.semilogy(isotherms, adiabat_p, color = 'gray', lw = 0.25, ls = '--')
    plt.semilogy(my_w, w_p, color = 'purple', lw = 0.5, ls = '--')
    w_labels = [w*1000. if w*1000 < 1 else int(w*1000) for w in target_w]
    for iw in range(len(target_w)):
        plt.text(my_w[np.where(w_p == 850)[0][0], iw], 850., str(w_labels[iw]), color = 'purple')
    plt.gca().invert_yaxis()
    plt.ylim([1050, 100])
    plt.gca().yaxis.grid()
    plt.yticks(wmo_p, wmo_p)
    plt.xlim([-40, 40])
#paper()

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
              CAPE = False, date = -999, my_title = '', temp_col = ['r'], 
              dew_col = ['b']):
    """
    Plots a skewT-logP diagram on the paper I've defined above.
    Requires an input temperature, dewpoint, pressure, and optional u- and v- 
    wind components. If CAPE is set to .True. then a call is made to calculate
    the surface-based parcel ascent. This provides CAPE, CIN computed using
    temperature.
    
    Here, CAPE is defined as the positively buoyant region between the level
    of free convection and the highest level of neutral buoyancy.
    
    CIN, however, is defined as the negatively buoyant region between the 
    lifting condensation level and the level of free convection.
    ---------------------------------------------------------------------------
    Optimal input units:
    temp    = deg C
    t_dew   = deg C
    p       = hPa
    u       = kts
    v       = kts
    """
    ### get the correct units for temp and t_dew ###
    # should be in deg C
    from scipy import interpolate 
    
    gs = gridspec.GridSpec(1, 2, width_ratios = [4, 1])
    gs.update(wspace=0.0, hspace=0.0)
    #paper()
    
    ## plot the skew T paper ##
    theta = dry_adiabats()
    thetae = moist_adiabats()
    my_w = const_ws()
    # skew the adiabats and the mixing ratio lines
    for sfc in range(len(T0)):
        theta[:,sfc] = skew(theta[:,sfc], adiabat_p)
    for sfc in range(len(T1)):
        thetae[:,sfc] = skew(thetae[:,sfc], adiabat_p)
    for sfc in range(len(target_w)):
        my_w[:,sfc] = skew(my_w[:,sfc], w_p)
    # get some isotherms
    therms = np.arange(-120., 40.1, 5.)
    isotherms = np.zeros((len(adiabat_p), len(therms)))
    for level in range(len(adiabat_p)):
        isotherms[level,:] = therms
    for sfc in range(len(therms)):
        isotherms[:,sfc] = skew(isotherms[:,sfc], adiabat_p)
    
    #start the plot for the paper
    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    ax = fig.add_subplot(gs[0])
    # dry adiabats
    ax.semilogy(theta, adiabat_p, color = 'g', lw = 0.25)
    # pseudo adiabats
    ax.semilogy(thetae, adiabat_p, color = 'b', lw = 0.25)
    # isotherms
    ax.semilogy(isotherms, adiabat_p, color = 'gray', lw = 0.25, ls = '--')
    # mixing ratios
    ax.semilogy(my_w, w_p, color = 'purple', lw = 0.5, ls = '--')
    # label the mixing ratios
    w_labels = [w*1000. if w*1000 < 1 else int(w*1000) for w in target_w]
    for iw in range(len(target_w)):
        ax.text(my_w[np.where(w_p == 850)[0][0], iw], 850., str(w_labels[iw]), color = 'purple')
    
    ## plot the sounding ##
    # make sure our input are lists
    if type(temp) != list:
        print 'WARNING! converting to a list of temperature profiles.'
        temp = [temp]
    if type(t_dew) != list:
        print 'WARNING! Converting to a list of dewpoint profiles.'
        t_dew = [t_dew]
    if type(p) != list:
        print 'WARNING! Converting to a list of pressure profiles.'
        p = [p]
    if type(u) != list:
        print 'WARNING! Converting to a list of uwind profiles.'
        u = [u]
    if type(v) != list:
        print 'WARNING! Converting to a list of vwind profiles.'
        v = [v]
    
    # some error checking to ensure that there is a matching temperature and dewpoint list
    if len(temp) != len(t_dew):
        print 'You have an unequal number of temperature and dewpoint profiles'
    if len(temp_col) != len(temp):
        print 'WARNING! The number of temperature colors provided does not match the number of temperature profiles. Taking the same dewpoint color for each profile.'
        while len(temp_col) != len(temp):
            temp_col.append(temp_col[-1])
    if len(dew_col) != len(t_dew):
        print 'WARNING! The number of dewpoint colors provided does not match the number of dewpoint profiles. Taking the same dewpoint color for each profile.'
        while len(dew_col) != len(t_dew):
            dew_col.append(dew_col[-1])
    
    for i in range(len(temp)):
        #each i in the list is a profile of t, td, and p
        # make sure they are the right units
        if (np.min(temp[i]) > 0) and (np.max(temp[i] > 100.)):
            temp[i] -= 273.15
        if (np.min(t_dew[i]) > 0) and (np.max(t_dew[i] > 100.)):
            t_dew[i] -= 273.15
        # plot the temperature
        print 'Plotting the temperature profile'
        #print 'Length temp[i] = ' + str(len(temp[i]))
        #print 'Length p[i] = ' + str(len(p[i]))
        ax.semilogy(skew(temp[i], p[i]), p[i], color = temp_col[i], lw = 2)
        # plot the dew point
        print 'Plotting the dewpoint profile'
        #print 'Length t_dew[i] = ' + str(len(temp[i]))
        ax.semilogy(skew(t_dew[i], p[i]), p[i], color = dew_col[i], lw = 2)
        if CAPE:
            my_q = getQ(t_dew[i]+273.15, np.zeros_like(t_dew[i])+100., p[i])
            my_CAPE, my_CIN, my_Parcel, LCLp, LFCp = getCAPE(temp[i], my_q, p[i], parcel_type = 1)
            my_title = 'CAPE = ' + str(round(my_CAPE, 0)) + 'Jkg$^{-1}$\n CIN = ' + str(round(my_CIN, 0)) + 'Jkg$^{-1}$'
            if np.nanmin(my_Parcel) > 0:
                my_Parcel = np.array([obs - 273.15 for obs in my_Parcel])
            ax.semilogy(skew(my_Parcel, p[i]), p[i], color = 'gray', lw = 2)
            ax.semilogy([-100, 100], [LCLp, LCLp], color = 'k')
            ax.semilogy([-100, 100], [LFCp, LFCp], color = 'k')
    ax.set_ylim([np.max([np.max([np.max(p_i) for p_i in p]), 1025.]), 100])
    ax.grid(axis = 'y')
    ax.set_yticks(wmo_p)
    ax.set_yticklabels(wmo_p)
    ax.set_xlim([-40, 40])
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Pressure (hPa)')
    
    if date != -999:
        my_title = date + '\n' + my_title
    ax.set_title(my_title)    
    
    ## plot the wind barbs ##
    for i in range(len(u)):
        if u[i][0] != -999:
            # windBarbs(u, v, p, gs)
            # only use wind obs on wmo pressure levels so we can see all the wind barbs
            #print 'Length p[i] = ' + str(len(p[i]))
            #print 'Length u[i] = ' + str(len(u[i]))
            u_wmo = interpolate.interp1d(p[i], u[i], fill_value = 'extrapolate')(wmo_p)
            v_wmo = interpolate.interp1d(p[i], v[i], fill_value = 'extrapolate')(wmo_p)
            
            # add wind barbs to the SkewT
            ax2 = fig.add_subplot(gs[1], sharey=ax)
            ax2.barbs(np.zeros_like(wmo_p), wmo_p, u_wmo, v_wmo, barbcolor = 'k', \
                clip_on = False, zorder=100)
            ax2.set_xlim([-1, 1])
            ax2.set_ylim([np.max([np.max([np.max(p_i) for p_i in p]), 1025.]), 100])
            print 'Plotted wind barbs'
            ax2.set_xticks([])
            ax2.axis('off')


    
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

