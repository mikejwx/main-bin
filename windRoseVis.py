import matplotlib as mpl
#Code to create wind roses for the station data and sounding data
#execfile('/home/xb899100/main/CloudTrails/CTID_vis.py')

#plt a wind rose for the cloud trail days
#reformat the testdates
MYDATES = {'Cloud Trails' : testDates,
           'Non-Trails' : testDates2,
           'Obscured' : testDates3,
           'Climatology' : testDates + testDates2 + testDates3}

for key in MYDATES.keys():
    myDates = []
    for date in MYDATES[key]:
        myDates.append(str(date))
    CL_SPD = []
    CL_DIR = []
    for eachDay in myDates:
        hour = eachDay.split()[1].split(':')[0]
        minute = eachDay.split()[1].split(':')[1]
        if minute == '15':
            minute = '20'
        elif minute == '45':
            minute = '50'
        #if hour == '12':
        eachDayindateTimeFormat = eachDay.split()[0] + ' ' + hour + ':' + minute + ':00'
        month = eachDay.split('-')[1]
        if (eachDayindateTimeFormat in data12['date-time']) and (month in ['06','07','08']):
            index = data12['date-time'].index(eachDayindateTimeFormat)
            CL_SPD.append(data12['spd'][index])
            CL_DIR.append(data12['DIR'][index])
            
    
    myres = 22.5 #This is commented out to use the myres in the CTID
    sectors = np.arange(0, 360, myres)
    speeds = np.arange(0, 41, 5)
    sorting = np.zeros((len(sectors), len(speeds)))
    calm = 0
    for ob in range(len(CL_SPD)):
        if CL_SPD[ob] > 0:
            i = int(np.round(CL_SPD[ob]/5, 0))
            myDir = CL_DIR[ob]
            if myDir > (360. - myres/2.):
                myDir -= 360.
            j = int(np.round(myDir/myres, 0))
            sorting[j, i] += 1
        elif CL_SPD == 0:
            calm += 1
            
    
    # Compute pie slices
    dirs = (sectors - myres/2.)*np.pi/180.
    fig = plt.gcf()
    fig.set_size_inches(8,8)
    ax = plt.subplot(111, projection='polar')
    cmap = mpl.cm.Reds
    norm = mpl.colors.BoundaryNorm(speeds, cmap.N)
    for spd in range(len(speeds)-1,0, -1):
        totals = np.zeros(len(sectors))
        for s in range(spd):
            totals += sorting[:,s]
        bars = ax.bar(dirs, totals/len(CL_SPD), width = myres*np.pi/180., bottom=0.0, color = cmap(speeds[spd]/40.))
    plt.ylim([0, 0.20])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    plt.title('Wind Rose for ' + key)
    plt.text(20*np.pi/180., 0.21, 'Mean Speed: ' + str(np.round(np.mean(CL_SPD),1)), fontsize = 18)
    ax1 = fig.add_axes([1, 0.15, 0.05, 0.8])
    cbar = mpl.colorbar.ColorbarBase(ax1, cmap = cmap, norm = norm, orientation = 'vertical')
    cbar.set_label('Wind Speed (ms $^{-1}$)')
    plt.savefig('/home/xb899100/main/WindRose_Surface_rwy12_'+key+'.png', bbox_inches = 'tight')
    plt.show()