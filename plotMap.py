import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

BDA_lon = -64.8
BDA_lat = 32.3
fig = plt.gcf()
fig.set_size_inches(12, 12)
m = Basemap(projection='mill',resolution='l', llcrnrlon=-180,
            llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90)
m.drawcoastlines()
m.fillcontinents(color='lightgreen',lake_color='cyan')
m.drawmapboundary(fill_color='aqua')
x,y = m(BDA_lon, BDA_lat)
m.plot(x,y, 'ro')


m = Basemap(projection='mill',resolution='h', llcrnrlon=-100,
            llcrnrlat=10,urcrnrlon=5,urcrnrlat=60)
fig = plt.gcf()
fig.set_size_inches(12, 12)
m.drawcoastlines()
m.fillcontinents(color='lightgreen',lake_color='cyan')
m.drawmapboundary(fill_color='aqua')
m.drawmeridians(np.arange(-180,180.,5), color = 'k', labels = [False, False, False, True])
m.drawparallels(np.arange(0.,90.,5), color = 'k', labels = [True, False, False, False])
x,y = m(BDA_lon, BDA_lat)
m.plot(x,y, 'ro')
lon = 0.
lon_min = -56.
lon_sec = -9.
lon = lon + (lon_min + (lon_sec/60.))/60.
lat = 51
lat_min = 26.
lat_sec = 22.
lat = lat + (lat_min + (lat_sec/60.))/60.
x, y = m(lon, lat)
m.plot(x, y, 'bo')

