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
data_dir = '/home/disk/r2d2/steed/cmip5/rcp8.5'

geo_em = '/home/disk/a125/steed/run/geo_em.d02.nc'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
vars_all = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
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
    
#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------
load_obs_all = 0
pickle_file = pickle_dir + '/obs_stats.pickle'

if load_obs_all == 0:
    #--- If not loading all obs, just load stats pickle file.
    if os.path.isfile(pickle_file):
        (obs,yy_obs) = pickle.load(open(pickle_file, 'rb'))
    else:
        sys.exit('cannot see ' + pickle_file)
else:
    #--- Load GHCND obs; get seasonal stats on them.
    elements = ['TMAX', 'TMIN', 'PRCP']
    (data_all, dts_all) = read_ghcnd(stns, elements, ghcnd_dir, sdt, edt)

    obs = {}
    yy_obs = {}
    
    var = 'TMAX'
    stat = 'max'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)
    var = 'TMAX'
    stat = 'avg'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)
    var = 'TMIN'
    stat = 'min'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)
    var = 'TMIN'
    stat = 'avg'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)
    var = 'PRCP'
    stat = 'tot'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)
    var = 'PRCP'
    stat = 'max'
    (obs[var,stat], yy_obs[var,stat]) = \
                    get_seasonal_stats_ghcnd(data_all,dts_all,stns, var,stat)

    pickle.dump((obs,yy_obs), open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
load_mod_all = 0
pickle_file = pickle_dir + '/mod_stats.pickle'

if load_mod_all == 0:
    #--- If not loading all model data, just load stats pickle file.
    if os.path.isfile(pickle_file):
        (mod,yy_mod) = pickle.load(open(pickle_file, 'rb'))
    else:
        sys.exit('cannot see ' + pickle_file)
else:
    pickle_file_all = pickle_dir + '/extract_' + sdt + '_' + edt + '.pickle'
    print 'reading pickle file ', pickle_file_all
    if os.path.isfile(pickle_file_all):
        (data_all, dts_unique, models, vars, stns) = \
                   pickle.load(open(pickle_file_all, 'rb'))
        print 'done!'
    else:
        (data_all, dts_unique, models, vars_all, stns) = \
                   load_extract_data(geo_em, stns, latpts, lonpts, elevs, \
                                     models, data_dir, sdt, edt, vars_all, \
                                     pickle_dir, lp)
    
    mod = {}
    yy_mod = {}

    var = 'PREC'
    stat = 'tot'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)
    var = 'PREC'
    stat = 'max'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)
    var = 'T2MAX'
    stat = 'max'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)
    var = 'T2MAX'
    stat = 'avg'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)
    var = 'T2MIN'
    stat = 'min'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)
    var = 'T2MIN'
    stat = 'avg'
    (mod[var,stat], yy_mod[var,stat]) = \
                    get_seasonal_stats(data_all,dts_unique,models,stns,var,stat)

    pickle.dump((mod,yy_mod), open(pickle_file,'wb'), -1)

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
    for s in range(len(seasons)):
        season = seasons[s]
#        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
#                    season + '_' + var + '_max.png'
#        titlein = 'Maximum 24-h ' + seasons_lab[season] + ' Precipitation, ' \
#                  + stn.upper() + ', ' + str(years[0]) + ' - ' + \
#                  str(years[len(years)-1])
#        iret = ts_mods(pcp_max, years, stn, models, titlein, plotfname, \
#                       var, 'max', mod_cols, [season])

        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                    season + '_' + var + '_tot.png'
        titlein = seasons_lab[season] + ' Total Precipitation, ' \
                  + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(pcp_tot, years, pcp_tot_obs, years_obs, stn, models, \
                       titlein, plotfname, var, 'tot', mod_cols, [season])
        sys.exit()
        
sys.exit()

##---------------------------------------------------------------------------
## Get seasonal temperature maxes and mins.
##---------------------------------------------------------------------------
#(tmx_max, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                      'T2MAX', 'max')
#(tmn_min, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                      'T2MIN', 'min')
#
##--- Do multiple stations on a plot for one model.
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

