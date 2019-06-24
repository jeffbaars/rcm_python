#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from igtg_check_utils import *
import matplotlib.pyplot as plt
#import pickle
from make_cmap import *

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir  = '/home/disk/spock/jbaars/rcm'
py_dir   = rcm_dir + '/python'
plot_dir = rcm_dir + '/plots'
#data_dir = rcm_dir + '/igtg.nml.works.3d.4d'
data_dir = '/home/disk/mass/steed/cmip5/rcp4.5/access1.0/common/' \
           'igtg.nml.works.3d.tos'
#data_dir = '/home/disk/mass/steed/cmip5/rcp4.5/access1.0/common/testdata'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#varname = 'ua'
varname = 'tos'
#file_old = data_dir + '/' + varname + '.1990010106.1991010100.nc'
#file_new = data_dir + '/' + varname + '_new.nc'
file_old = data_dir + '/' + 'indata'
file_new = data_dir + '/' + 'outdata'

itime = 1
ilev  = 1

title_old = varname + ' OLD, itime = ' + str(itime) + ', ilev = ' + \
            str(ilev)
plotfname_old = './' + varname + '_old_itime' + str(itime) + '_ilev' + \
                str(ilev) + '.png'
title_new = varname + ' NEW, itime = ' + str(itime) + ', ilev = ' + \
            str(ilev)
plotfname_new = './' + varname + '_new_itime' + str(itime) + '_ilev' + \
                str(ilev) + '.png'

#minlat = 20.0
#maxlat = 50.0
#minlon = -120.0
#maxlon = -50.0

minlat = 20.0
maxlat = 70.0
minlon = -175.0
maxlon = -85.0

#---------------------------------------------------------------------------
# Make a color map using make_cmap.
#---------------------------------------------------------------------------
if (varname == 'tos'):
    levs_c = levs_tos
    cmap_c = cmap_tos
elif(varname == 'ua'):
    levs_c = levs_spd
    cmap_c = cmap_spd
        
levs_norm_c = []
for i in range(len(levs_c)):
    x = float(levs_c[i])
    norm_c = (x - min(levs_c)) / (max(levs_c) - min(levs_c))
    levs_norm_c.append(norm_c)

my_cmap = make_cmap(cmap_c, bit = True, position = levs_norm_c)

#---------------------------------------------------------------------------
# Load old and new files.
#---------------------------------------------------------------------------
(lat_old, lon_old, data_old) = load_nc(file_old, varname, itime, ilev)

print 'lat_old, lon_old, data_old = ', lat_old, lon_old, data_old

(lat_new, lon_new, data_new) = load_nc(file_new, varname, itime, ilev)

print 'lat_new, lon_new, data_new = ', lat_new, lon_new, data_new

#sys.exit()

        
#---------------------------------------------------------------------------
# Plot old data.
#---------------------------------------------------------------------------
mapper(lat_old, lon_old, data_old, minlat, maxlat, minlon, maxlon, \
       levs_c, my_cmap, maplw, colorbar_lab, title_old, plotfname_old)

mapper(lat_new, lon_new, data_new, minlat, maxlat, minlon, maxlon, \
       levs_c, my_cmap, maplw, colorbar_lab, title_new, plotfname_new)

sys.exit()

