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

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
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
#vars_mod =['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
vars_mod = ['PREC', 'T2MAX', 'T2MIN']

sdt = '197001'
edt = '209912'

var_obs_dict = {
    'PREC': 'PRCP',
    'T2MAX': 'TMAX',
    'T2MIN': 'TMIN'
    }
stats_dict = {
    'PREC': ['max', 'tot'],
    'T2MAX': ['max', 'avg'],
    'T2MIN': ['min', 'max', 'avg']
    }

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']

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

#---------------------------------------------------------------------------
# Load GHCND obs.
#---------------------------------------------------------------------------
load_obs_all = 0
pickle_file = pickle_dir + '/obs_stats.pickle'
if load_obs_all == 0:
    #--- If not loading all obs, just load stats pickle file.
    if os.path.isfile(pickle_file):
        print 'reading pickle file ', pickle_file
        (obs,yy_obs) = pickle.load(open(pickle_file, 'rb'))
    else:
        sys.exit('cannot see ' + pickle_file)
else:
    vars_obs = []
    for var in vars_mod:
        vars_obs.append(var_obs_dict[var])

    #--- Load GHCND obs; get seasonal stats on them.
    (data_all, dts_all) = read_ghcnd(stns, vars_obs, ghcnd_dir, sdt, edt)
    obs = {}
    yy_obs = {}
    for var in vars_mod:
        varo = var_obs_dict[var]
        stats = stats_dict[var]
        for stat in stats:
            (obs[varo,stat], yy_obs[varo,stat]) = get_seasonal_stats_ghcnd(\
                data_all,dts_all,stns,varo,stat)

    pickle.dump((obs,yy_obs), open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
# Load model data.
#---------------------------------------------------------------------------
load_mod_all = 0
pf_mod_stats = pickle_dir + '/mod_stats.pickle'
if load_mod_all == 0:
    #--- If not loading all model data, just load stats pickle file.
    if os.path.isfile(pf_mod_stats):
        print 'reading pickle file ', pf_mod_stats
        (mod,yy_mod) = pickle.load(open(pf_mod_stats, 'rb'))
    else:
        sys.exit('cannot see ' + pf_mod_stats)
else:

    for m in range(len(models)):
        model = models[m]
        pf = pickle_dir + '/extract_' + model +'_' + sdt + '_' + edt + '.pickle'
        if not os.path.isfile(pf):
            (data_all, dts_unique, models, vars_all, stns) = \
                       load_extract_data(geo_em, stns, latpts, lonpts, elevs, \
                                         models, data_dirs, sdt, edt, \
                                         vars_mod, pickle_dir, 0)

    data_all = {}
    dts_all = []
    for m in range(len(models)):
        model = models[m]
        pf = pickle_dir + '/extract_' + model +'_' + sdt + '_' + edt + '.pickle'
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
                data_all,dts_unique,models,stns,var,stat)
    pickle.dump((mod,yy_mod), open(pf_mod_stats,'wb'), -1)

#sys.exit()

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
#    sys.exit()
sys.exit()

