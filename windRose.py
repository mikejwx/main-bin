import matplotlib as mpl
#Code to create wind roses for the station data and sounding data
execfile('/home/xb899100/bin/projectFun.py') #base functions for my project
#execfile('/home/xb899100/bin/readSoundingData.py')
#execfile('/home/xb899100/main/CloudTrails/CTID.py')

def getWindRoseData(dates, months = ['05', '06', '07', '08', '09', '10']):
    #reformat the testdates
    #month can be a list of month number strings e.g. ['05', '06'] for may and june
    myDates = []
    for date in dates:
        myDates.append(str(date).split()[0])
    
    CL_SPD = []
    CL_DIR = []
    for eachDay in range(len(soundingdate)):
        day = str(soundingdate[eachDay]).split()[0]
        time = str(soundingdate[eachDay]).split()[1]
        month = day.split('-')[1]
        if (day in myDates) and (time == '12:00:00') and (month in months):
            U, DIR = fromComponents(interpolated['u_p'][:,eachDay], interpolated['v_p'][:,eachDay], True)
            Umean, DIRmean = presWeight(950., 850., P, U, DIR, 'm/s')    
            CL_SPD.append(Umean)
            CL_DIR.append(DIRmean)
    
    #myres = 10. #This is commented out to use the myres in the CTID
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
    return sectors, speeds, sorting, dirs, CL_SPD, CL_DIR
#==============================================================================
# #fig = plt.gcf()
# #fig.set_size_inches(8,8)
# #ax = plt.subplot(111, projection='polar')
# cmap = mpl.cm.Reds
# norm = mpl.colors.BoundaryNorm(speeds, cmap.N)
# for spd in range(len(speeds)-1,0, -1):
#     totals = np.zeros(len(sectors))
#     for s in range(spd):
#         totals += sorting[:,s]
#     bars = ax.bar(dirs, totals/len(CL_SPD), width = myres*np.pi/180., bottom=0.0, color = cmap(speeds[spd]/40.))
# plt.ylim([0, 0.33])
# ax.set_theta_zero_location("N")
# ax.set_theta_direction(-1)
# #plt.title('Wind Rose for JJA Cloud Trails')
# plt.text(45*np.pi/180., 0.38, 'Calm: ' + str(calm), fontsize = 18)
# ax1 = fig.add_axes([1, 0.15, 0.05, 0.8])
# cbar = mpl.colorbar.ColorbarBase(ax1, cmap = cmap, norm = norm, orientation = 'vertical')
# cbar.set_label('Wind Speed (ms $^{-1}$)')
#==============================================================================
