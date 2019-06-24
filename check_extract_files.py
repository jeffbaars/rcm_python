#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
from utils_ghcnd_obs import *
import pickle
import csv
import json

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
#plot_dir   = rcm_dir + '/plots'
plot_dir   = rcm_dir + '/plots/new'
pickle_dir = rcm_dir + '/pickle'
data_dir   = rcm_dir + '/data'
obs_dir    = rcm_dir + '/obs'
ghcnd_dir  = obs_dir + '/ghcnd_all'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

geo_em = data_dir + '/geo_em.d02.nc'
station_file = obs_dir + '/station_file.txt'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#vars_mod = ['PREC', 'T2MAX', 'T2MIN']
vars_mod = ['T2MAX']

sdt = '197001'
edt = '209912'

#--- GHCND data downloaded through June 2019, so 2019 is a partial year, so
#--- only include obs through 2018.
edt_obs = '201812' 

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

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']
models = ['access1.0']

temp_max = 273 + 60
temp_min = 273 - 60

for m in range(len(models)):
    data_all = {}
    dts_all = []
    model = models[m]

    #--- Get list of extracted files for this model.
    print 'model = ', model
    files = get_extract_files(data_dirs, model, sdt, edt, 0)
    files = list(sorted(set(files)))

    #--- Loop over all extract files, reading them in.
    for f in range(len(files)):
        (vardat, ntimes, dt) = load_wrfout(files[f], vars_mod)

        for v in range(len(vars_mod)):
            var = vars_mod[v]
            dat_c = vardat[var]
            x,y,z = dat_c.shape
            for k in range(z):
                max_c = np.amax(dat_c[:,:,k])
                if max_c < temp_min or max_c > temp_max:
                    print model, ', ', files[f], ', ', k
                          
                





'''
#---------------------------------------------------------------------------
# Load model data.
#---------------------------------------------------------------------------
load_mod_all = 1
pf_mod_stats = pickle_dir + '/mod_stats.pkl'
if load_mod_all == 0:
    #--- If not loading all model data, just load stats pickle file.
    if os.path.isfile(pf_mod_stats):
        print 'reading pickle file ', pf_mod_stats
        (mod,yy_mod) = pickle.load(open(pf_mod_stats, 'rb'))
    else:
        sys.exit('cannot see ' + pf_mod_stats)
else:

    (data_all, dts_unique, models, vars_all, stns) = \
               load_extract_data(geo_em, stns, latpts, lonpts, elevs, models, \
                                 data_dirs, sdt, edt, vars_mod, pickle_dir, 0)
    sys.exit()
    
    data_all = {}
    dts_all = []
    for m in range(len(models)):
        model = models[m]
        pf = get_mod_pkl_name(pickle_dir, model, sdt, edt)
        print 'loading ', pf
        (model, data_c, dts_unique, vars, stns) = pickle.load(open(pf, 'rb'))
        data_all[model] = data_c
        dts_all = dts_all + dts_unique
    mod = {}
    yy_mod = {}
    for var in vars_mod:
        stats = stats_dict[var]        
        for stat in stats:
            (mod[var,stat], yy_mod[var,stat]) = get_seasonal_stats(\
                data_all, dts_unique, models, stns, var, stat)
    pickle.dump((mod,yy_mod), open(pf_mod_stats,'wb'), -1)

sys.exit()

#---------------------------------------------------------------------------
# Make a time series plots of each station, variable, stat and season.
#---------------------------------------------------------------------------
#--- Do multiple models on a plot for one station.
for ss in range(len(stns)):
    stn = stns[ss]
    for v in range(len(vars_mod)):
        var = vars_mod[v]
        varo = var_obs_dict[var]
        stats = stats_dict[var]
        for st in range(len(stats)):
            stat = stats[st]
            for s in range(len(seasons)):
                season = seasons[s]
                plotfname = get_plotfname(plot_dir,stn,sdt,edt,season,var,stat)
                titlein = get_title(season,stn,var,stat,yy_mod[var,stat])
                iret = ts_mods(mod[var,stat],yy_mod[var,stat], \
                               obs[varo,stat],yy_obs[varo,stat], stn, \
                               models, titlein, plotfname, var, stat, \
                               mod_cols, [season])
#                sys.exit()
    sys.exit()
sys.exit()

'''
