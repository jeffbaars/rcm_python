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
    
##---------------------------------------------------------------------------
## Load GHCND obs; get seasonal stats on them.
##---------------------------------------------------------------------------
#elements = ['TMAX', 'TMIN', 'PRCP']
#
#load_pickle = 1
#pickle_file = pickle_dir + '/ghcnd_test.pickle'
#if load_pickle and os.path.isfile(pickle_file):
#    print 'reading pickle file ', pickle_file
#    (data_all, dts_all, elements, stns) = \
#               pickle.load(open(pickle_file, 'rb'))
#    print 'done!'
#else:
#    (data_all, dts_all) = read_ghcnd(stns, elements, ghcnd_dir, sdt, edt)
#    pickle.dump((data_all, dts_all, elements, stns), \
#                open(pickle_file,'wb'), -1)
#
#var = 'TMAX'
#(tmx_mx_obs, yyyy_obs) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                                  stns, var, 'max')
#(tmx_avg, years) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                            stns, var, 'avg')
#var = 'TMIN'
#(tmn_mn, years) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                           stns, var, 'min')
#(tmn_avg, years) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                            stns, var, 'avg')
#var = 'PRCP'
#(pcp_tot_obs, years_obs) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                                    stns, var, 'tot')
#(pcp_max_obs, years_obs) = get_seasonal_stats_ghcnd(data_all, dts_all, \
#                                                    stns, var, 'max')

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
led = 0
if led == 1:
    lp = 1
    (data_all, dts_unique, models, vars_all, stns) = \
               load_extract_data(geo_em, stns, latpts, lonpts, elevs, models, \
                                 data_dir, sdt, edt, vars_all, pickle_dir, lp)

pickle_file = pickle_dir + '/pcp_tot.pickle'
if os.path.isfile(pickle_file):
    (pcp_tot, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (pcp_tot, years) = get_seasonal_stats(data_all, dts_unique, models, \
                                          stns, 'PREC', 'tot')
    pickle.dump((pcp_tot,years), open(pickle_file,'wb'), -1)

pickle_file = pickle_dir + '/pcp_max.pickle'
if os.path.isfile(pickle_file):
    (pcp_max, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (pcp_max, years) = get_seasonal_stats(data_all, dts_unique, models, \
                                          stns, 'PREC', 'max')
    pickle.dump((pcp_max,years), open(pickle_file,'wb'), -1)
    
pickle_file = pickle_dir + '/tmx_max.pickle'
if os.path.isfile(pickle_file):
    (tmx_max, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (tmx_max, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                          'T2MAX', 'max')
    pickle.dump((tmx_max,years), open(pickle_file,'wb'), -1)
    
pickle_file = pickle_dir + '/tmn_min.pickle'
if os.path.isfile(pickle_file):
    (tmn_min, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (tmn_min, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                          'T2MIN', 'min')
    pickle.dump((tmn_min,years), open(pickle_file,'wb'), -1)

pickle_file = pickle_dir + '/tmx_avg.pickle'
if os.path.isfile(pickle_file):
    (tmx_avg, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (tmx_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                          'T2MAX', 'avg')
    pickle.dump((tmx_avg,years), open(pickle_file,'wb'), -1)


pickle_file = pickle_dir + '/tmn_avg.pickle'
if os.path.isfile(pickle_file):
    (tmn_avg, years) = pickle.load(open(pickle_file, 'rb'))
else:
    (tmn_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                          'T2MIN', 'avg')
    pickle.dump((tmn_avg,years), open(pickle_file,'wb'), -1)

sys.exit()

#---------------------------------------------------------------------------
# Get seasonal precipitation totals and make a time series plot of them.
#---------------------------------------------------------------------------
var = 'PREC'
(pcp_tot, years) = get_seasonal_stats(data_all, dts_unique, models, \
                                      stns, var, 'tot')
(pcp_max, years) = get_seasonal_stats(data_all, dts_unique, models, \
                                      stns, var, 'max')

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
        titlein = 'Total 24-h ' + seasons_lab[season] + ' Precipitation, ' \
                  + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(pcp_tot, years, pcp_tot_obs, years_obs, stn, models, \
                       titlein, plotfname, var, 'tot', mod_cols, [season])
        sys.exit()
        
sys.exit()

#---------------------------------------------------------------------------
# Get seasonal temperature maxes and mins.
#---------------------------------------------------------------------------
(tmx_max, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MAX', 'max')
(tmn_min, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MIN', 'min')

#--- Do multiple stations on a plot for one model.
mod = 'miroc5'
for s in range(len(seasons)):
    season = seasons[s]

    #--- Max temperature.
    var = 'T2MAX'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_max.png'
    titlein = 'Maximum Seasonal Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(tmx_max, years, stns, mod, titlein, plotfname, var, 'max', \
              stn_cols, [season])

    #--- Min temperature.
    var = 'T2MIN'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_min.png'
    titlein = 'Minimum Seasonal Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(tmn_min, years, stns, mod, titlein, plotfname, var, 'min', \
              stn_cols, [season])

sys.exit()

#--- Do multiple models on a plot for one station.
for ss in range(len(stns)):
    stn = stns[ss]
    for s in range(len(seasons)):
        season = seasons[s]

        #--- Max temperature.
        var = 'T2MAX'
        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                    season + '_' + var + '_max.png'
        titlein = 'Maximum Seasonal Temperature, ' + seasons_lab[season] + \
                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(tmx_max, years, stn, models, titlein, plotfname, \
                       var, 'max', stn_cols, [season])

        #--- Min temperature.
        var = 'T2MIN'
        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                    season + '_' + var + '_min.png'
        titlein = 'Minimum Seasonal Temperature, ' + seasons_lab[season] + \
                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(tmn_min, years, stn, models, titlein, plotfname, \
                       var, 'min', stn_cols, [season])

#---------------------------------------------------------------------------
# Get seasonal temperature averages.
#---------------------------------------------------------------------------
(tmx_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MAX', 'avg')
(tmn_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MIN', 'avg')

#--- Do multiple stations on a plot for one model.
mod = 'gfdl-cm3'
for s in range(len(seasons)):
    season = seasons[s]

    #--- Max temperature.
    var = 'T2MAX'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_avg.png'
    titlein = 'Average Maximum Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(tmx_avg, years, stns, mod, titlein, plotfname, var, 'avg', \
              stn_cols, [season])

    #--- Min temperature.
    var = 'T2MIN'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_avg.png'
    titlein = 'Average Minimum Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(tmn_avg, years, stns, mod, titlein, plotfname, var, 'avg', \
              stn_cols, [season])


#--- Do multiple models on a plot for one station.
for ss in range(len(stns)):
    stn = stns[ss]
    for s in range(len(seasons)):
        season = seasons[s]

        #--- Max temperature.
        var = 'T2MAX'
        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                    season + '_' + var + '_avg.png'
        titlein = 'Average Maximum Temperature, ' + seasons_lab[season] + \
                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(tmx_avg, years, stn, models, titlein, plotfname, \
                       var, 'avg', stn_cols, [season])

        #--- Min temperature.
        var = 'T2MIN'
        plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                    season + '_' + var + '_avg.png'
        titlein = 'Average Minimum Temperature, ' + seasons_lab[season] + \
                  ', ' + stn.upper() + ', ' + str(years[0]) + ' - ' + \
                  str(years[len(years)-1])
        iret = ts_mods(tmn_avg, years, stn, models, titlein, plotfname, \
                       var, 'avg', stn_cols, [season])

sys.exit()

#---------------------------------------------------------------------------
# Wind Speed.
#---------------------------------------------------------------------------
mod = 'gfdl-cm3'
var = 'SPDUV10MEAN'
(spdavg_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                         var, 'avg')
for s in range(len(seasons)):
    season = seasons[s]

    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_avg.png'
    titlein = 'Average Wind Speed, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(spdavg_avg, years, stns, mod, titlein, plotfname, var, \
                   'avg', stn_cols, [season])

var = 'SPDUV10MAX'
(spdmax_max, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                         var, 'max')
for s in range(len(seasons)):
    season = seasons[s]

    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '_max.png'
    titlein = 'Maximum Wind Speed, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts_stns(spdmax_max, years, stns, mod, titlein, plotfname, var, \
                   'max', stn_cols, [season])


sys.exit()

