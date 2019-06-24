#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import pickle
from utils_date import *
import csv

ll_lat =   45.3
ur_lat =   49.3
ll_lon = -125.0
ur_lon = -116.5

fs      = 9
titlefs = 9
width   = 10
height  = 8
#width   = 30
#height  = 24
maplw   = 1.0

geo_em = '/home/disk/a125/steed/run/geo_em.d02.nc'
mask_file = '/home/disk/spock/jbaars/rcm/data/mask_pts.dat'

#---------------------------------------------------------------------------
# Read mask lat/lon points file.
#---------------------------------------------------------------------------
mask = open(mask_file)
csvReader = csv.reader(mask)
maski = []
maskj = []
for row in csvReader:
    maski.append(int(row[0]))
    maskj.append(int(row[1]))

#---------------------------------------------------------------------------
# Load geo_em file.
#---------------------------------------------------------------------------
print 'Loading ', geo_em
nc = NetCDFFile(geo_em, 'r')
lat = nc.variables['XLAT_M'][0,:,:]
lon = nc.variables['XLONG_M'][0,:,:]
hgt = nc.variables['HGT_M'][0,:,:]    
lat = np.transpose(lat)
lon = np.transpose(lon)
hgt = np.transpose(hgt)
nc.close()

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

res = 'i'
#res = 'h'
        
fig = plt.figure(figsize=(width,height))
# left, bottom, width, height:
ax = fig.add_axes([0.00,0.05,0.99,0.91])
map = Basemap(resolution = res,projection='lcc',\
              llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
              urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
              lat_0=lat_ctr,lon_0=lon_ctr,lat_1=(ur_lat - ll_lat))

#--- Get lat and lon data in map's x/y coordinates.
x,y = map(lon, lat)

#--- Draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth = maplw)
map.drawstates(linewidth = maplw)
map.drawcountries(linewidth = maplw)

#--- Draw the edge of the map projection region (the projection limb)
map.drawmapboundary(linewidth = maplw)

#-- Add dots grid points.
(nx,ny) = lon.shape

#for x in range(nx):
#    for y in range(ny):
#        lat_c = lat[x,y]
#        lon_c = lon[x,y]
#        if (lat_c >= ll_lat and lat_c <= ur_lat and \
#            lon_c >= ll_lon and lon_c <= ur_lon):
#            xs, ys = map(lon_c, lat_c)
#            plt.plot(xs, ys, 'ro', markersize = 3)
#            plt.text(xs, ys, str(x) + ',' + str(y), fontsize=8)

plt.title('WA state mask grid points', fontsize=titlefs, fontweight='bold')

for m in range(len(maski)):
    lat_c = lat[maski[m],maskj[m]]
    lon_c = lon[maski[m],maskj[m]]
    xs, ys = map(lon_c, lat_c)
    plt.plot(xs, ys, 'ro', markersize = 4)


#--- Save plot.
plotfname = './WA_mask.png'
print 'xli ', plotfname, ' &'
plt.savefig(plotfname)

plt.close()

sys.exit()

