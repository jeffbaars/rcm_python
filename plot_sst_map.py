#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from collections import defaultdict
from utils_cmap import *
import pickle
from utils_date import *
from utils_wrfout import *
from matplotlib import colors as c
import numpy.ma as ma

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
model = 'ccsm4'
dt = '1999072900'
#geo_em = '/home/disk/spock/jbaars/rcm/data/geo_em.d02.nc'
geo_em = '/home/disk/riker/steed/geo_em.d02.nc.urban.water.nc'

wrfout = './wrf2d.d02.1hr.' + dt + '.nc'

gcm = '/home/disk/r2d2/steed/cmip5/rcp8.5/ccsm4/data/ncungrib/' + \
      'ncFILE:1999-07-01_00'

##----------------------------------------------------------------------------
## Read gcm.
##----------------------------------------------------------------------------
#nc = NetCDFFile(gcm, 'r')
#lat_gcm = nc.variables['lat'][:]
#lon_gcm = nc.variables['lon'][:]
#lm_gcm  = nc.variables['LANDSEA'][:,:]
#nc.close()
#
#xs, ys = np.meshgrid(lon_gcm, lat_gcm)
#
##----------------------------------------------------------------------------
## Plot it.
##----------------------------------------------------------------------------
#titlein = 'GCM Land Mask, ' + model + ', ' + dt
#plotfname = './gcm_landmask.png'
#zoom = 'Z1'
#
#ur_lat = maxlat[zoom]
#ur_lon = maxlon[zoom]
#ll_lat = minlat[zoom]
#ll_lon = minlon[zoom]
#lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
#lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)
#
#if (zoom == 'Z1'):
#    res = 'i'
#elif (zoom == 'Z2'):
#    res = 'h'
#        
#fig = plt.figure(figsize=(width,height))
## left, bottom, width, height:
#ax = fig.add_axes([0.00,0.05,0.99,0.91])
#map = Basemap(resolution = res,projection='lcc',\
#              llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
#              urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
#              lat_0=lat_ctr,lon_0=lon_ctr,lat_1=(ur_lat - ll_lat))
#
##--- Get lat and lon data in map's x/y coordinates.
##x,y = map(lon_gcm, lat_gcm)
#x,y = map(xs, ys)
#
##--- Draw coastlines, country boundaries, fill continents.
#map.drawcoastlines(linewidth = maplw)
#map.drawstates(linewidth = maplw)
#map.drawcountries(linewidth = maplw)
#
##--- Draw the edge of the map projection region (the projection limb)
#map.drawmapboundary(linewidth = maplw)
##--- Draw lat/lon grid lines every 30 degrees.
#map.drawmeridians(np.arange(0, 360, 30), linewidth = maplw)
#map.drawparallels(np.arange(-90, 90, 30), linewidth = maplw)
#
#cMap = c.ListedColormap(['cornflowerblue','darkkhaki'])
#plt.pcolormesh(x, y, lm_gcm, cmap=cMap)
#
#plt.title(titlein, fontsize=titlefs, fontweight='bold')
#
##--- Save plot.
#print 'xli ', plotfname, ' &'
#plt.savefig(plotfname)
#
#plt.close()
#
#sys.exit()

#----------------------------------------------------------------------------
# Read geo_em.
#----------------------------------------------------------------------------
nc = NetCDFFile(geo_em, 'r')
lat  = nc.variables['XLAT_M'][0,:,:]
lon  = nc.variables['XLONG_M'][0,:,:]
hgt  = nc.variables['HGT_M'][0,:,:]
lm   = nc.variables['LANDMASK'][0,:,:]
luse = nc.variables['LU_INDEX'][0,:,:]

lat  = np.transpose(lat)
lon  = np.transpose(lon)
hgt  = np.transpose(hgt)
lm   = np.transpose(lm)
luse = np.transpose(luse)
nc.close()

urban = 13
water = 17

luse[luse == urban] = 100
luse[luse == water] = 101
luse[luse < 100]    = 102

#----------------------------------------------------------------------------
# Read wrfout.
#----------------------------------------------------------------------------
nc = NetCDFFile(wrfout, 'r')
sst = nc.variables['SST'][0,:,:]
sst = np.transpose(sst)
nc.close()

sst = ((sst - 273.15) * 9.0/5.0) + 32

#----------------------------------------------------------------------------
# Plot it.
#----------------------------------------------------------------------------
plotfname = './test.png'
titlein = 'Surface Temperature, Land Mask, Urban Land Use, ' + model + \
          ', ' + dt
zoom = 'Z2'

minlatZ2 = 46.5
maxlatZ2 = 49.0
minlonZ2 = -124.0
maxlonZ2 = -121.0
#minlatZ2 = 44.0
#maxlatZ2 = 47.0
#minlonZ2 = -125.0
#maxlonZ2 = -120.0

#ur_lat = maxlat[zoom]
#ur_lon = maxlon[zoom]
#ll_lat = minlat[zoom]
#ll_lon = minlon[zoom]
ur_lat = maxlatZ2
ur_lon = maxlonZ2
ll_lat = minlatZ2
ll_lon = minlonZ2
lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

if (zoom == 'Z1'):
    res = 'i'
elif (zoom == 'Z2'):
    res = 'h'
        
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
#--- Draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 30), linewidth = maplw)
map.drawparallels(np.arange(-90, 90, 30), linewidth = maplw)

cMap = c.ListedColormap(['darkgray', 'royalblue','darkkhaki'])
CS = plt.contour(x, y, sst, [52,56,60,64,68,72,76,80], color = 'white')
plt.clabel(CS, inline=1, fontsize=10, fmt='%2.0f')
plt.pcolormesh(x, y, luse, cmap=cMap)

#--- Plot station locations.
stn_plot = ['KPDX', 'KSEA']
stn_lats = [45.59083, 47.44472]
stn_lons = [-122.60028, -122.31361]
for sp in range(len(stn_plot)):
    xs, ys = map(stn_lons[sp], stn_lats[sp])
    plt.plot(xs, ys, 'rs', markersize = 5, markeredgewidth = 1.1)
    plt.text(xs, ys, stn_plot[sp])

plt.title(titlein, fontsize=titlefs, fontweight='bold')

#--- Save plot.
print 'xli ', plotfname, ' &'
plt.savefig(plotfname)
#plt.show()

plt.close()
