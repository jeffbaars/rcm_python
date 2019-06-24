#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
#from make_cmap import *
#from netCDF4 import Dataset as NetCDFFile
#from utils_wrfout import *
#from utils_date import *
#import pickle

#from ghcndextractor import *
#from ghcndextractor import ghcndextractor
import codecs
from utils_ghcnd_obs import *

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
#ghcnFolder  = rcm_dir + '/obs'
ghcnd_dir   = rcm_dir + '/obs/ghcnd_all'
data_dir   = '/home/disk/r2d2/steed/cmip5/rcp8.5'

geo_em = '/home/disk/a125/steed/run/geo_em.d02.nc'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#sdt = '202001'
#edt = '202912'
sdt = '197001'
edt = '209912'

dtfmt = "%Y%m%d%H"

stns   = ['KSEA', 'KSMP', 'KYKM', 'KGEG', 'KPDX', 'KMFR', 'KHQM']
models = ['gfdl-cm3', 'miroc5']
stn_cols = ['b', 'r', 'g', 'c', 'm', 'k', 'grey']
mod_cols = ['b', 'r']
latpts = [47.44472, 47.27667, 46.56417, 47.62139, 45.59083, \
          42.38111, 46.97278]
lonpts = [-122.31361, -121.33722, -120.53361, -117.52778, -122.60028, \
          -122.87222, -123.93028]
elevs  = [130.0, 1207.0, 333.0, 735.0, 8.0, 405.0, 4.0]

ghcnd_dict = {
    'KSEA': 'USW00024233',
    'KSMP': 'USW00024237',
    'KYKM': 'USW00024243',
    'KGEG': 'USW00024157',
    'KPDX': 'USW00024229',
    'KMFR': 'USW00024225',
    'KHQM': 'USW00094225'}

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
elements = ['TMAX', 'TMIN', 'PRCP']

for s in range(len(stns)):
    stn_ghcnd = ghcnd_dict[stns[s]]
    file_c = ghcnd_dir + '/' + stn_ghcnd + '.dly'

    if os.path.isfile(file_c):
        print 'reading ', file_c
        (datag, dts_all) = read_ghcnd(file_c, elements)

    test = get_seasonal_stats_ghcnd(datag, dts_all)

#    for dt in dts_all:
#        print dt
#
#    for key in datag.keys():
#        print(key)
#
#    sys.exit()
#
#    el = 'TMAX'
#    dt = '19701205'
#    key = (el,dt)
##    if (el,yr,mo,dy) in datag:
#    if key in datag:        
#        print 'datag[key] = ', datag[key]
#    else:
#        print 'cannot see key in datag! key = ', key
#    sys.exit()

sys.exit()


