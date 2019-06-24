#!/usr/bin/python
import  os, os.path, time, glob, re, math
import sys, string, readline, h5py
from Scientific.IO.NetCDF import *
import numpy as np
#from datetime import datetime, timedelta, time
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
#from rw_plotsettings import *
import pickle
import csv
from numpy import genfromtxt

maplw = 0.8

plotfname = './test_pts2.png'
title = 'pts'
#ptsfile = '/home/disk/mass/jbaars/cmip5/rcp8.5/miroc5/test_unstaggertos/pts.txt'
ptsfile = '/home/disk/mass/jbaars/cmip5/rcp8.5/miroc5/test_unstaggertos/pts3.txt'

indices = []
lats = []
lons = []
f = open(ptsfile, 'rt')
try:
    reader = csv.reader(f)
    for index, lat, lon in reader:
        indices.append(float(index))
        lats.append(float(lat))
        lons.append(float(lon))
finally:
    f.close()

maxlat = max(lats) + 0.5
minlat = min(lats) - 0.5
maxlon = max(lons) + 0.5
minlon = min(lons) - 0.5

ur_lat = maxlat
ur_lon = maxlon
ll_lat = minlat
ll_lon = minlon
lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

width   = 12
height  = 4
fs      = 9
titlefs = 9

fig = plt.figure(figsize=(width,height))
# left, bottom, width, height:
ax = fig.add_axes([0.00,0.05,0.99,0.91])

#---
map = Basemap(resolution='i',projection='lcc',\
              llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
              urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
              lat_0=lat_ctr,lon_0=lon_ctr,lat_1=(ur_lat - ll_lat))

#--- Draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth = maplw)
map.drawstates(linewidth = maplw)
map.drawcountries(linewidth = maplw)

#--- Draw the edge of the map projection region (the projection limb)
map.drawmapboundary(linewidth = maplw)

#--- Draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0, 360, 30), linewidth = maplw)
map.drawparallels(np.arange(-90, 90, 30), linewidth = maplw)

#-- Add dots for grid points for 2nd grid.
for i in range(len(lats)):

    xs, ys = map(lons[i], lats[i])

    if(indices[i] == 2):
        mark = 'gx'
        ms = 5.0
        plt.plot(xs, ys, mark, markersize = 3, mew=ms)
        xs2 = xs
        ys2 = ys
    if (indices[i] == 1):
        mark = 'ro'
        ms = 1.0
        plt.plot([xs,xs2], [ys,ys2], linewidth=2, color='g')
    else:
        mark = 'bo'
        ms = 1.0

    plt.plot(xs, ys, mark, markersize = 6, mew=ms)

#    strp = "%d,%d,%.3f" % (ii+1, jj+1, tos_new[ii+1,jj+1])
#    plt.text(xs, ys, strp, fontsize=6)

plt.title(title, fontsize=titlefs, fontweight='bold')

#--- Save plot.
print 'creating ', plotfname
plt.savefig(plotfname)

plt.close()

