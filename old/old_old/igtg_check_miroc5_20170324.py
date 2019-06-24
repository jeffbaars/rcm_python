#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from igtg_check_utils_miroc5 import *
import matplotlib.pyplot as plt
#import pickle
from make_cmap import *
from Scientific.IO.NetCDF import *

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir  = '/home/disk/spock/jbaars/rcm'
py_dir   = rcm_dir + '/python'
plot_dir = rcm_dir + '/plots'
data_dir = '/home/disk/mass/jbaars/cmip5/rcp8.5/miroc5/test_unstaggertos'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
varname = 'tos'
file_old = data_dir + '/' + 'inlatlon'
file_new = data_dir + '/' + 'indata'

itime = 1
ilev  = 1

title = 'both'
plotfname = './both3.png'
#title = 'inlatlon'
#plotfname = './inlatlon.png'
#title = 'indata'
#plotfname = './indata.png'

i1 = 4
j1 = 107
i2 = 174
j2 = 183

#minlat = 50.0
#maxlat = 70.0
#minlon = -10.0
#maxlon = 20.0
#minlat = 55.0
#maxlat = 65.0
#minlon = 0.0
#maxlon = 10.0

#---------------------------------------------------------------------------
# Load old and new files.
#---------------------------------------------------------------------------
(lat_old, lon_old) = load_nc(file_old, 0)
(lat_new, lon_new) = load_nc(file_new, 1)

print 'lat_old = ', lat_old
print 'lat_new = ', lat_new
print ''
print 'lon_old = ', lon_old  
print 'lon_new = ', lon_new  

#sys.exit()

#nx,ny = np.shape(lon_old)
#print 'nx, ny = ', nx, ny
#for ii in range(nx):
#    for jj in range(ny):
#
#        if (lat_old[ii,jj] > minlat and lat_old[ii,jj] < maxlat):
#
#            lon_c = lon_old[ii,jj]
#            if (lon_c > 180.0):
#                lon_c = lon_c - 360
#            
#            print 'is ', minlon, ' < ', lon_c, ' < ', maxlon
#            if (lon_c > minlon and lon_c < maxlon):
#                print 'yes'
#            else:
#                print 'no'
#sys.exit()


lat_thresh = 2.0
lon_thresh = 3.0
print 'lat_old[i1,j1], lon_old[i1,j1] = ', lat_old[i1,j1], lon_old[i1,j1]
maxlat = lat_old[i1,j1] + lat_thresh
minlat = lat_old[i1,j1] - lat_thresh
maxlon = lon_old[i1,j1] + lon_thresh
minlon = lon_old[i1,j1] - lon_thresh

#---------------------------------------------------------------------------
# Plot old data.
#---------------------------------------------------------------------------
mapper(i1, j1, i2, j2, lat_old, lon_old, lat_new, lon_new, \
       minlat, maxlat, minlon, maxlon, title, plotfname)

sys.exit()

