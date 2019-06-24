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
#levs_spd = [ 0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5,
#             25.0, 27.5, 30.0, 32.5, 35.0, 37.5 ]
levs_spd = [ -24, -21, -18, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 18, \
             21, 24]
cmap_spd = [(255,255,255), \
            (192,192,192), \
            (128,128,128 ), \
            (0,255,255), \
            (32,178,170), \
            (0,255,0), \
            (0,128,0), \
            (255,0,204), \
            (199,21,133), \
            (0,0,255), \
            (0,0,128), \
            (255,255,0), \
            (255,204,17), \
            (255,69,0), \
            (0,0,0),\
            (255,255,255)]

#---------------------------------------------------------------------------
# Load geog file.
#---------------------------------------------------------------------------
def load_nc(filein, varname, itime, ilev):
    file = NetCDFFile(filein, 'r')
    lat_tmp  = file.variables['lat'][:]
    lon_tmp  = file.variables['lon'][:]
    data_all = file.variables[varname]
    fillval = data_all._FillValue[0]
    data_tmp = file.variables[varname][itime,ilev,:,:]
    file.close()

    #--- Convert from 0-to-360 longitudes to -180-to-180.
    ipos = np.where(lon_tmp >= 180.0)
    iwest = ipos[0]
    lon_tmp[iwest] = lon_tmp[iwest] - 360.0

    #--- transpose data array to be x/y (lon/lat).
    dataout = np.transpose(data_tmp)
    dataout[dataout == fillval] = np.NAN

    #--- Fill in 2-D lat / lon terr arrays.
    nx = len(lon_tmp)
    ny = len(lat_tmp)

    lonout = np.empty([nx,ny])
    latout = np.empty([nx,ny])
    for y in range(ny):
        latout[:,y] = lat_tmp[y]
    for x in range(nx):
        lonout[x,:] = lon_tmp[x]

    return latout, lonout, dataout

#---------------------------------------------------------------------------
# Mapper.
#---------------------------------------------------------------------------
def mapper(lat, lon, grid, minlat, maxlat, minlon, maxlon, \
           levs, cmap, maplw, colorbar_lab, title, plotfname):

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
    map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
#    map = Basemap(resolution='i',projection='lcc',\
#                  llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
#                  urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
#                  lat_0=lat_ctr,lon_0=lon_ctr,lat_1=(ur_lat - ll_lat))

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
        
    cs = plt.contourf(x, y, grid, levs, cmap=cmap)

    cbar = map.colorbar(cs, location='bottom', pad="3%", size=0.1, ticks=levs)
    cbar.set_label(colorbar_lab, fontsize=fs, size=fs-1)
    cbar.ax.tick_params(labelsize=fs-1)

    plt.title(title, fontsize=titlefs, fontweight='bold')

    #--- Save plot.
    print 'creating ', plotfname
    plt.savefig(plotfname)

    plt.close()

