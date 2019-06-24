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

colorbar_lab = 'm/s';
maplw = 0.8

#---------------------------------------------------------------------------
# Load 2-D lat/lon netcdf data file.
#---------------------------------------------------------------------------
def load_nc(filein, latlon2d):

    print 'reading ', filein
    f = NetCDFFile(filein, 'r')
    ndims = len(f.dimensions)
    dims  = f.dimensions
    dimnames = f.dimensions.keys()

    for d in range(ndims):
        print d, dimnames[d], dims[dimnames[d]]
        #--- Find ny dimension.
        if (dimnames[d] == 'lat' or dimnames[d] == 'rlat'):
            ny = dims[dimnames[d]]
        #--- Find nx dimension.
        if (dimnames[d] == 'lon' or dimnames[d] == 'rlon'):
            nx = dims[dimnames[d]]

    if (latlon2d == 1):
        lon_tmp = f.variables['lon'][:,:]
        lat_tmp = f.variables['lat'][:,:]
        lonout = np.transpose(lon_tmp)
        latout = np.transpose(lat_tmp)
    else:
        lon_tmp = f.variables['lon'][:]
        lat_tmp = f.variables['lat'][:]
        lonout = np.empty([nx,ny])
        latout = np.empty([nx,ny])
        for y in range(ny):
            latout[:,y] = lat_tmp[y]
        for x in range(nx):
            lonout[x,:] = lon_tmp[x]

    f.close()

    return latout, lonout

#---------------------------------------------------------------------------
# Mapper.
#---------------------------------------------------------------------------
def mapper(i1, j1, i2, j2, lat1, lon1, lat2, lon2,
           minlat, maxlat, minlon, maxlon, title, plotfname):

    ur_lat = maxlat
    ur_lon = maxlon
    ll_lat = minlat
    ll_lon = minlon
    lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
    lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

    width   = 6.5
    height  = 6.5
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

#    #--- Get lat and lon data in map's x/y coordinates.
#    x1,y1 = map(lon1, lat1)
#    x2,y2 = map(lon2, lat2)    

    #-- Add dots for grid points for 1st grid.
    nx,ny = np.shape(lon1)
    print 'nx, ny = ', nx, ny
    for ii in range(nx):
        for jj in range(ny):
            if (lat1[ii,jj] > minlat and lat1[ii,jj] < maxlat):
                lon_c = lon1[ii,jj]
                if (lon_c > 180.0):
                    lon_c = lon_c - 360
                if (lon_c > minlon and lon_c < maxlon):
                    if (ii == i1 and jj == j1):
                        mark = 'ro'
                    else:
                        mark = 'bo'                        
                    xs, ys = map(lon1[ii,jj], lat1[ii,jj])
                    plt.plot(xs, ys, mark, markersize = 5)
                    plt.text(xs, ys, str(ii) + ',' + str(jj), fontsize=7)
                    
    #-- Add dots for grid points for 2nd grid.
    nx,ny = np.shape(lon2)
    print 'nx, ny = ', nx, ny
    for ii in range(nx):
        for jj in range(ny):
            if (lat2[ii,jj] > minlat and lat2[ii,jj] < maxlat):
                lon_c = lon2[ii,jj]
                if (lon_c > 180.0):
                    lon_c = lon_c - 360
                if (lon_c > minlon and lon_c < maxlon):
                    if (ii == i2 and jj == j2):
                        mark = 'rx'
                        ms = 2.0
                    else:
                        mark = 'kx'
                        ms = 1.2
                    xs, ys = map(lon2[ii,jj], lat2[ii,jj])
                    plt.plot(xs, ys, mark, markersize = 6, mew=ms)
                    plt.text(xs, ys, str(ii) + ',' + str(jj), fontsize=6)

    plt.title(title, fontsize=titlefs, fontweight='bold')

    #--- Save plot.
    print 'creating ', plotfname
    plt.savefig(plotfname)

    plt.close()

