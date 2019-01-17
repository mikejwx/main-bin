#Plot the Bermuda DEM

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

#read in the DEM
DEM = Dataset('/home/xb899100/bin/BermudaDEM.nc')
DEMlon = DEM.variables['longitude'][:]
DEMlat = DEM.variables['latitude'][:]
dem = DEM.variables['elevation'][0,:,:]

#set up the map space
myMap = Basemap(projection = 'mill', resolution = 'h', llcrnrlon = min(DEMlon),
                llcrnrlat = min(DEMlat), urcrnrlon = max(DEMlon), 
                urcrnrlat = max(DEMlat))
x,y = np.meshgrid(DEMlon, DEMlat)
X,Y = myMap(x, y)

#set up levels for the color bar
topoLevels = np.linspace(0., 100., 101)

#plot a 9x6" figure with colors centred on 0, lat/lon every 0.1 degrees and 
#key ticks every 20 m
plt.clf()
fig = plt.gcf()
fig.set_size_inches(9,6)
myMap.contourf(X, Y, dem, levels = topoLevels, cmap = plt.cm.gist_earth, 
               vmin = -np.max(dem), vmax = np.max(dem))
myMap.drawparallels(np.arange(30.,34.,0.1), color = 'cyan', 
                    labels = [True, False, False, False])
myMap.drawmeridians(np.arange(-75.,360.,0.1), color = 'cyan', 
                    labels = [False, False, False, True])
plt.colorbar(ticks = [0., 20., 40., 60., 80., 100.], label = 'elevation (m)')
plt.savefig('BermudaDEM.png')
