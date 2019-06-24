#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
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
vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#sdt = '202001'
#edt = '202912'
sdt = '197001'
edt = '209912'
models = ['gfdl-cm3']
dtfmt = "%Y%m%d%H"

stns   = ['KSEA', 'KSMP', 'KYKM']
stn_cols = ['b', 'r', 'g']
latpts = [47.44472, 47.27667, 46.56417]
lonpts = [-122.31361, -121.33722, -120.53361]
elevs  = [130.0, 1207.0, 333.0]
    
##---------------------------------------------------------------------------
## Test.
##---------------------------------------------------------------------------
#x = range(0,len(vars),2)
#print 'x = ', x
#
##test = vars[x]
##print 'test = ', test
#
#test = [ vars[i] for i in x ]
#print 'test = ', test
#
#sys.exit()

#ndays = {}
#season = 'winter'
#yyyy = '1970'
#ndays[season+yyyy] = 26.0
#if (season+yyyy in ndays):
#    print 'hi'
#else:
#    print 'bye'
#sys.exit()

#season = 'winter'
#seasonsin = ['winter', 'spring']
#if (season in seasonsin):
#    print 'hi'
#else:
#    print 'bye'
#
#sys.exit()

#---------------------------------------------------------------------------
# Loop over each models, grabbing data from extract wrfout files.
#---------------------------------------------------------------------------
load_pickle = 1
pickle_file = pickle_dir + '/extract_' + sdt + '_' + edt + '.pickle'
if load_pickle and os.path.isfile(pickle_file):
    (data_all, dts_unique, models, vars, stns) = \
               pickle.load(open(pickle_file, 'rb'))
else:

    #--- Load geo em file to get lat/lon/elev grids.  Also interpolate model
    #--- elevation to stns lat/lon point locations.
    print 'Loading:\n', geo_em
    (lat, lon, hgt) = load_geo_em(geo_em)
    elevs_mod  = []
    elevs_diff = []    
    for ns in range(len(stns)):
        dat_interp = bilinear_interpolate(lat, lon, hgt, \
                                          latpts[ns], lonpts[ns], 0)
        elevs_mod.append(dat_interp)
        elevs_diff.append(elevs[ns] - dat_interp)

    data_all = {}
    dts_all = []
    for m in range(len(models)):

        model = models[m]
        data_dir_c = data_dir + '/' + model + '/extract'

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dir_c, model, sdt, edt)

        print 'Loading: '
        #--- Loop over all extract files, reading them in.
        for f in range(len(files)):

            print files[f]
            (vardat, ntimes, dt) = load_wrfout(files[f], vars)

            #--- Loop over each vars, interpolating data to stns for
            #--- each ntimes.
            for v in range(len(vars)):

                var = vars[v]
                vardat_c = vardat[var]
                
                #--- Loop over each ntime, grabbing interpolated data for
                #--- variable var for each stns.
                for n in range(ntimes):
                    dt_c = time_increment(dt, n, dtfmt)
                    for ns in range(len(stns)):
                        dat_interp = bilinear_interpolate(lat, lon, \
                                                          vardat_c[:,:,n],\
                                                          latpts[ns], \
                                                          lonpts[ns], 0)

                        #--- Do lapse rate correction for temperature.
                        if (var == 'T2MAX' or var == 'T2MIN'):
                            dat_c = dat_interp - std_lapse * elevs_diff[ns]
                        else:
                            dat_c = dat_interp

                        data_all[model, var, stns[ns], dt_c] = dat_c
                        dts_all.append(dt_c)

    #--- Get sorted unique list of all dates seen.
    dts_unique = list(sorted(set(dts_all)))
    
    #--- Create pickle file.
    pickle.dump((data_all, dts_unique, models, vars, stns), \
                open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
# Get seasonal precipitation totals and make a time series plot of them.
#---------------------------------------------------------------------------
var = 'PREC'
(pcp_tot, years) = get_seasonal_stats(data_all, dts_unique, models, \
                                      stns, var, 'tot')

mod = 'gfdl-cm3'
for s in range(len(seasons)):
    season = seasons[s]
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '.png'
    titlein = 'Total ' + seasons_lab[season] + ' Precipitation, ' + \
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts(pcp_tot, years, stns, mod, titlein, plotfname, 'PREC', stn_cols, \
              [season])

#---------------------------------------------------------------------------
# Get seasonal temperature averages.
#---------------------------------------------------------------------------
(tmx_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MAX', 'avg')
(tmn_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                      'T2MIN', 'avg')

mod = 'gfdl-cm3'
for s in range(len(seasons)):
    season = seasons[s]

    #--- Max temperature.
    var = 'T2MAX'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '.png'
    titlein = 'Average Maximum Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts(tmx_avg, years, stns, mod, titlein, plotfname, var, \
              stn_cols, [season])

    #--- Min temperature.
    var = 'T2MIN'
    plotfname = plot_dir + '/' + mod + '_' + sdt + '_' + edt + '_' + season + \
                '_' + var + '.png'
    titlein = 'Average Minimum Temperature, ' + seasons_lab[season] + ', ' +\
              mod.upper() + ', ' + str(years[0]) + ' - ' + \
              str(years[len(years)-1])
    iret = ts(tmn_avg, years, stns, mod, titlein, plotfname, var, \
              stn_cols, [season])

#---------------------------------------------------------------------------
# Wind Speed.
#---------------------------------------------------------------------------
(spdavg_avg, years) = get_seasonal_stats(data_all, dts_unique, models, stns, \
                                         'SPDUV10MEAN', 'avg')

sys.exit()

