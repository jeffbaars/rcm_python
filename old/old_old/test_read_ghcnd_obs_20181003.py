#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
#from make_cmap import *
#from netCDF4 import Dataset as NetCDFFile
#from utils_wrfout import *
#from utils_date import *
import pickle

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

#dtfmt = "%Y%m%d%H"

stns   = ['KSEA', 'KSMP', 'KYKM', 'KGEG', 'KPDX', 'KMFR', 'KHQM']
#stns   = ['KSEA']
models = ['gfdl-cm3', 'miroc5']
stn_cols = ['b', 'r', 'g', 'c', 'm', 'k', 'grey']
mod_cols = ['b', 'r']
latpts = [47.44472, 47.27667, 46.56417, 47.62139, 45.59083, \
          42.38111, 46.97278]
lonpts = [-122.31361, -121.33722, -120.53361, -117.52778, -122.60028, \
          -122.87222, -123.93028]
elevs  = [130.0, 1207.0, 333.0, 735.0, 8.0, 405.0, 4.0]

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
elements = ['TMAX', 'TMIN', 'PRCP']

load_pickle = 1
pickle_file = pickle_dir + '/ghcnd_test.pickle'
if load_pickle and os.path.isfile(pickle_file):
    print 'reading pickle file ', pickle_file
    (data_all, dts_all, elements, stns) = \
               pickle.load(open(pickle_file, 'rb'))
    print 'done!'
else:
    (data_all, dts_all) = read_ghcnd(stns, elements, ghcnd_dir, sdt, edt)
    pickle.dump((data_all, dts_all, elements, stns), \
                open(pickle_file,'wb'), -1)

var = 'TMAX'
(tmx_mx, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'max')
(tmx_avg, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'avg')
var = 'TMIN'
(tmn_mn, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'min')
(tmn_avg, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'avg')
var = 'PRCP'
(pcp_mx, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'max')
(pcp_tot, years) = get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, 'tot')

for year in years:
    key = ('fall','KSEA',year)
    if (key in pcp_tot):
        print year, pcp_tot[key]
    else:
        print key, ' not in there'

sys.exit()


