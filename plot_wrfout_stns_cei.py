#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from utils_wrfout import *
from utils_date import *
from utils_ghcnd_obs import *
from utils_cei import *
import cPickle as pickle  #-- cPickle is faster, but deprecated in python 3.0
import csv

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
#plot_dir   = rcm_dir + '/plots'
plot_dir   = rcm_dir + '/plots/new'
pickle_dir = rcm_dir + '/pickle'
data_dir   = rcm_dir + '/data'
ghcnd_dir  = data_dir + '/ghcnd/ghcnd_all'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

geo_em = data_dir + '/geo_em.d02.nc'
station_file = data_dir + '/station_file.txt'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#vars_mod =['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
vars_mod = ['PREC', 'T2MAX', 'T2MIN', 'T2MEAN']

sdt = '197001'
#edt = '209912'
edt = '197912'

#--- GHCND data downloaded through June 2019, so 2019 is a partial year, so
#--- only include obs through 2018.
#edt_obs = '201812'
edt_obs = '197912' 

var_obs_dict = {
    'PREC': 'PRCP',
    'T2MAX': 'TMAX',
    'T2MIN': 'TMIN',
    'T2MEAN': 'T2MEAN'    
    }
stats_dict = {
    'PREC': ['max', 'tot'],
    'T2MAX': ['max', 'avg'],
    'T2MIN': ['min', 'max', 'avg'],
    'T2MEAN': ['min', 'max', 'avg'],    
    }

#models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
#          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
#          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']
models = ['access1.0', 'access1.3']

mod_cols = ['indigo', 'blue', 'deepskyblue', \
            'darkgreen', 'lime', 'yellow', \
            'magenta', 'red', 'salmon', 'gray', 'darkgray', 'lightblue']

#---------------------------------------------------------------------------
# Load stations file.
#---------------------------------------------------------------------------
stns   = []
latpts = []
lonpts = []
elevs  = []
with open(station_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        stns.append(row[0])
        latpts.append(float(row[1]))
        lonpts.append(float(row[2]))
        elevs.append(float(row[3]))

##---------------------------------------------------------------------------
## Load GHCND obs.
##---------------------------------------------------------------------------
#vars_obs = []
#for var in vars_mod:
#    vars_obs.append(var_obs_dict[var])
#
#pf = pickle_dir + '/obs_daily_test.pkl'
#if os.path.isfile(pf):
#    print 'Loading ', pf
#    (obs_all, obs_dts_all, stns, latpts, lonpts, elevs) = \
#               pickle.load(open(pf, 'rb'))
#else:
#    (obs_all, obs_dts_all) = read_ghcnd(stns, vars_obs, ghcnd_dir, sdt, edt_obs)
#    print 'Creating ', pf
#    pickle.dump((obs_all, obs_dts_all, stns, latpts, lonpts, elevs), \
#                open(pf,'wb'), -1)
#
#---------------------------------------------------------------------------
# Load model data.
#---------------------------------------------------------------------------
#pf_all = pickle_dir + '/mod_daily.pkl'
pf_all = 'junk'
if os.path.isfile(pf_all):
    print 'Loading ', pf_all
    (mod_all, mod_dts_all) = pickle.load(open(pf_all, 'rb'))
else:
    mod_all = {}
    mod_dts_all = []
    for m in range(len(models)):
        model = models[m]
        pf = get_mod_pkl_name(pickle_dir, model, sdt, edt)
        print 'Loading ', pf
        (model, data_c, dts_unique, vars, stns) = pickle.load(open(pf, 'rb'))
        mod_all[model] = data_c
        mod_dts_all = mod_dts_all + dts_unique
#    print 'Creating ', pf_all
#    pickle.dump((mod_all, mod_dts_all), open(pf_all,'wb'), -1)

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
pf = pickle_dir + '/cei_daily_perc_mxmn.pkl'
if os.path.isfile(pf):
    (daily_perc_mxmn) = pickle.load(open(pf, 'rb'))
else:
    daily_perc_mxmn = get_daily_perc_mxmn(models, stns, mod_all, sdt)
    print 'Creating ', pf
    pickle.dump((daily_perc_mxmn), open(pf,'wb'), -1)

test = get_monthly_mxmn(daily_perc_mxmn, mod_all, models, stns, sdt)

sys.exit()
