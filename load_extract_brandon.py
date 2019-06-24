#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
#from utils_ghcnd_obs import *
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
vars_mod =['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#vars_mod = ['PREC', 'T2MAX', 'T2MIN']

sdt = '197001'
edt = '197003'
#sdt = '197001'
#edt = '209912'

models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', \
            'access1.3', 'canesm2', 'noresm1-m', \
            'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2' ]
mod_cols = ['indigo', 'blue', 'deepskyblue', \
            'darkgreen', 'lime', 'yellow', \
            'magenta', 'red', 'salmon']
stn_cols = ['b', 'r', 'g', 'c', 'm', 'k', 'grey']

pf_tmp = pickle_dir + '/ext_MODEL_' + sdt + '_' + edt + '.pickle'

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
# Create pickle file of data for each model.
#---------------------------------------------------------------------------
(data_all, dts_unique, models, vars_all, stns) = \
           load_extract_data(geo_em, stns, latpts, lonpts, elevs, \
                             models, data_dirs, sdt, edt, \
                             vars_mod, pickle_dir, 0)

#---------------------------------------------------------------------------
# Load pickle file for each model.
#---------------------------------------------------------------------------
data_all = {}
dts_all = []
for m in range(len(models)):
    model = models[m]
    pf = pf_tmp.replace('MODEL', model )       
    print 'loading ', pf
    (model, data_c, dts_unique, vars, stns) = pickle.load(open(pf, 'rb'))
    data_all[model] = data_c
    dts_all = dts_all + dts_unique

##---------------------------------------------------------------------------
## Make a time series plots of each station, variable, stat and season.
##---------------------------------------------------------------------------
##--- Do multiple models on a plot for one station.
#for ss in range(len(stns)):
#    stn = stns[ss]
#    for v in range(len(vars_mod)):
#        var = vars_mod[v]
#        varo = var_obs_dict[var]
#        stats = stats_dict[var]
#        for st in range(len(stats)):
#            stat = stats[st]
#            for s in range(len(seasons)):
#                season = seasons[s]
#                plotfname = get_plotfname(plot_dir,stn,sdt,edt,season,var,stat)
#                titlein = get_title(season,stn,var,stat,yy_mod[var,stat])
#                iret = ts_mods(mod[var,stat],yy_mod[var,stat], \
#                               obs[varo,stat],yy_obs[varo,stat], stn, \
#                               models, titlein, plotfname, var, stat, \
#                               mod_cols, [season])

sys.exit()

