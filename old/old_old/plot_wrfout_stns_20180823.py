#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
from utils_ghcnd_obs import *
import pickle

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
ghcnd_dir  = rcm_dir + '/obs/ghcnd_all'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

geo_em = '/home/disk/a125/steed/run/geo_em.d02.nc'

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

models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', 'access1.3']
mod_cols = ['b',        'r',      'g',          'y']

stns   = ['KSEA', 'KSMP', 'KYKM', 'KGEG', 'KPDX', 'KMFR', 'KHQM']
stn_cols = ['b', 'r', 'g', 'c', 'm', 'k', 'grey']
latpts = [47.44472, 47.27667, 46.56417, 47.62139, 45.59083, \
          42.38111, 46.97278]
lonpts = [-122.31361, -121.33722, -120.53361, -117.52778, -122.60028, \
          -122.87222, -123.93028]
elevs  = [130.0, 1207.0, 333.0, 735.0, 8.0, 405.0, 4.0]
    
#---------------------------------------------------------------------------
# Load obs.
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
    #--- Load GHCND obs; get seasonal stats on them.
    (data_all, dts_all) = read_ghcnd(stns, vars_obs, ghcnd_dir, sdt, edt)
    obs = {}
    yy_obs = {}
    for var in vars_obs:
        for stat in stats:
            (obs[var,stat], yy_obs[var,stat]) = get_seasonal_stats_ghcnd(\
                data_all,dts_all,stns,var,stat)

    pickle.dump((obs,yy_obs), open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
# Load model data.
#---------------------------------------------------------------------------
load_mod_all = 0
pickle_file = pickle_dir + '/mod_stats.pickle'

if load_mod_all == 0:
    #--- If not loading all model data, just load stats pickle file.
    if os.path.isfile(pickle_file):
        print 'reading pickle file ', pickle_file
        (mod,yy_mod) = pickle.load(open(pickle_file, 'rb'))
    else:
        sys.exit('cannot see ' + pickle_file)
else:
    pickle_file_all = pickle_dir + '/extract_' + sdt + '_' + edt + '.pickle'
    pickle_file_all = 'testjunk'
    print 'reading pickle file ', pickle_file_all
    if os.path.isfile(pickle_file_all):
        (data_all, dts_unique, models, vars, stns) = \
                   pickle.load(open(pickle_file_all, 'rb'))
        print 'done!'
    else:
        (data_all, dts_unique, models, vars_all, stns) = \
                   load_extract_data(geo_em, stns, latpts, lonpts, elevs, \
                                     models, data_dirs, sdt, edt, vars_mod, \
                                     pickle_dir, 0)
#        sys.exit()
        
    mod = {}
    yy_mod = {}
    for var in vars_mod:
        stats = stats_dict[var]        
        for stat in stats:
            (mod[var,stat], yy_mod[var,stat]) = get_seasonal_stats(\
                data_all,dts_unique,models,stns,var,stat)

    pickle.dump((mod,yy_mod), open(pickle_file,'wb'), -1)

#sys.exit()
#
#var = 'PREC'
#stat = 'tot'
#key = (var,stat)
#test = mod[key]
#years = yy_mod[key]
#
#season = 'winter'
##model = 'bcc-csm1.1'
#model = 'gfdl-cm3'
#stn = 'KSEA'
#
#for yyyy in years:
#    key = (season,model,stn,yyyy)
#    dat_c = test[key]
#    print yyyy, dat_c


#for key in mod:
#    print key, mod[key]

#sys.exit()

#---------------------------------------------------------------------------
# Get seasonal precipitation totals and make a time series plot of them.
#---------------------------------------------------------------------------
##--- Do multiple stations on a plot for one model.
#mod = 'miroc5'
#for s in range(len(seasons)):
#    season = seasons[s]
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_tot.png'
#    titlein = 'Total ' + seasons_lab[season] + ' Precipitation, ' + \
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(pcp_tot, years, stns, mod, titlein, plotfname, 'PREC', \
#                   'tot', stn_cols, [season])
#
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_max.png'
#    titlein = 'Maximum 24-h ' + seasons_lab[season] + ' Precipitation, ' + \
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(pcp_max, years, stns, mod, titlein, plotfname, 'PREC', \
#                   'max', stn_cols, [season])

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
#            sys.exit()
sys.exit()

#---------------------------------------------------------------------------
# Get seasonal temperature maxes and mins.
#---------------------------------------------------------------------------
#--- Do multiple stations on a plot for one model.
#mod = 'miroc5'
#for s in range(len(seasons)):
#    season = seasons[s]
#
#    #--- Max temperature.
#    var = 'T2MAX'
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_max.png'
#    titlein = 'Maximum Seasonal Temperature, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(tmx_max, years, stns, mod, titlein, plotfname, var, 'max', \
#              stn_cols, [season])
#
#    #--- Min temperature.
#    var = 'T2MIN'
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_min.png'
#    titlein = 'Minimum Seasonal Temperature, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(tmn_min, years, stns, mod, titlein, plotfname, var, 'min', \
#              stn_cols, [season])
#
#sys.exit()
#
##--- Do multiple models on a plot for one station.
#for ss in range(len(stns)):
#    stn = stns[ss]
#    for s in range(len(seasons)):
#        season = seasons[s]
#
#        #--- Max temperature.
#        var = 'T2MAX'
#        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
#                    season + '_' + var + '_max.png'
#        titlein = 'Maximum Seasonal Temperature, ' + seasons_lab[season] + \
#                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
#                  str(years[len(years)-1])
#        iret = ts_mods(tmx_max, years, stn, models, titlein, plotfname, \
#                       var, 'max', stn_cols, [season])
#
#        #--- Min temperature.
#        var = 'T2MIN'
#        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
#                    season + '_' + var + '_min.png'
#        titlein = 'Minimum Seasonal Temperature, ' + seasons_lab[season] + \
#                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
#                  str(years[len(years)-1])
#        iret = ts_mods(tmn_min, years, stn, models, titlein, plotfname, \
#                       var, 'min', stn_cols, [season])
#
##---------------------------------------------------------------------------
## Get seasonal temperature averages.
##---------------------------------------------------------------------------
#(tmx_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                      'T2MAX', 'avg')
#(tmn_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                      'T2MIN', 'avg')
#
##--- Do multiple stations on a plot for one model.
#mod = 'gfdl-cm3'
#for s in range(len(seasons)):
#    season = seasons[s]
#
#    #--- Max temperature.
#    var = 'T2MAX'
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_avg.png'
#    titlein = 'Average Maximum Temperature, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(tmx_avg, years, stns, mod, titlein, plotfname, var, 'avg', \
#              stn_cols, [season])
#
#    #--- Min temperature.
#    var = 'T2MIN'
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_avg.png'
#    titlein = 'Average Minimum Temperature, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(tmn_avg, years, stns, mod, titlein, plotfname, var, 'avg', \
#              stn_cols, [season])
#
#
##--- Do multiple models on a plot for one station.
#for ss in range(len(stns)):
#    stn = stns[ss]
#    for s in range(len(seasons)):
#        season = seasons[s]
#
#        #--- Max temperature.
#        var = 'T2MAX'
#        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
#                    season + '_' + var + '_avg.png'
#        titlein = 'Average Maximum Temperature, ' + seasons_lab[season] + \
#                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
#                  str(years[len(years)-1])
#        iret = ts_mods(tmx_avg, years, stn, models, titlein, plotfname, \
#                       var, 'avg', stn_cols, [season])
#
#        #--- Min temperature.
#        var = 'T2MIN'
#        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
#                    season + '_' + var + '_avg.png'
#        titlein = 'Average Minimum Temperature, ' + seasons_lab[season] + \
#                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
#                  str(years[len(years)-1])
#        iret = ts_mods(tmn_avg, years, stn, models, titlein, plotfname, \
#                       var, 'avg', stn_cols, [season])
#
#sys.exit()
#
##---------------------------------------------------------------------------
## Wind Speed.
##---------------------------------------------------------------------------
#mod = 'gfdl-cm3'
#var = 'SPDUV10MEAN'
#(spdavg_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                         var, 'avg')
#for s in range(len(seasons)):
#    season = seasons[s]
#
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_avg.png'
#    titlein = 'Average Wind Speed, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(spdavg_avg, years, stns, mod, titlein, plotfname, var, \
#                   'avg', stn_cols, [season])
#
#var = 'SPDUV10MAX'
#(spdmax_max, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                         var, 'max')
#for s in range(len(seasons)):
#    season = seasons[s]
#
#    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
#                '_' + var + '_max.png'
#    titlein = 'Maximum Wind Speed, ' + seasons_lab[season] + ', ' +\
#              mod.upper() + ', ' + str(years[0]) + ' - ' + \
#              str(years[len(years)-1])
#    iret = ts_stns(spdmax_max, years, stns, mod, titlein, plotfname, var, \
#                   'max', stn_cols, [season])
#
#
#sys.exit()

