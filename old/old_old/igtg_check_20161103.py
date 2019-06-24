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
data_dir = rcm_dir + '/igtg.nml.works.3d.4d'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
varname = 'ua'
file_old = data_dir + '/' + varname + '.1990010106.1991010100.nc'
file_new = data_dir + '/' + varname + '_new.nc'

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

minlat = 20.0
maxlat = 50.0
minlon = -120.0
maxlon = -50.0

#minlat = 20.0
#maxlat = 90.0
#minlon = -175.0
#maxlon = -85.0

#---------------------------------------------------------------------------
# Make a color map using make_cmap.
#---------------------------------------------------------------------------
levs_norm_spd = []
for i in range(len(levs_spd)):
    x = float(levs_spd[i])
    norm_c = (x - min(levs_spd)) / (max(levs_spd) - min(levs_spd))
    levs_norm_spd.append(norm_c)

my_cmap = make_cmap(cmap_spd, bit = True, position = levs_norm_spd)

#---------------------------------------------------------------------------
# Load old and new files.
#---------------------------------------------------------------------------
(lat_old, lon_old, data_old) = load_nc(file_old, varname, itime, ilev)
(lat_new, lon_new, data_new) = load_nc(file_new, varname, itime, ilev)
        
#---------------------------------------------------------------------------
# Plot old data.
#---------------------------------------------------------------------------
#mapper(lat_old, lon_old, data_old, minlat, maxlat, minlon, maxlon, \
#       levs_spd, my_cmap, maplw, colorbar_lab, title_old, plotfname_old)
#
#sys.exit()

mapper(lat_new, lon_new, data_new, minlat, maxlat, minlon, maxlon, \
       levs_spd, my_cmap, maplw, colorbar_lab, title_new, plotfname_new)

sys.exit()

