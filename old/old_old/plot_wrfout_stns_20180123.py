#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
#import matplotlib.pyplot as plt
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
sdt = '202001'
edt = '202912'
models = ['gfdl-cm3']
dtfmt = "%Y%m%d%H"

stns = ['KSEA', 'KSMP', 'KYKM']
latpts = [47.44472, 47.27667, 46.56417]
lonpts = [-122.31361, -121.33722, -120.53361]

#---------------------------------------------------------------------------
# Loop over each models, grabbing data from extract wrfout files.
#---------------------------------------------------------------------------
pickle_file = pickle_dir + '/extract_' + sdt + '_' + edt + '.pickle'
if os.path.isfile(pickle_file):
    (data_all, dts_unique, models, vars, stns) = \
               pickle.load(open(pickle_file, 'rb'))
else:
    print 'Loading ', geo_em
    (lat, lon) = load_geo_em(geo_em)

    data_all = {}
    dts_all = []
    for m in range(len(models)):

        model = models[m]
        data_dir_c = data_dir + '/' + model + '/extract'

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dir_c, model, sdt, edt)

        #--- Loop over all extract files, reading them in.
        for f in range(len(files)):

            print 'Loading ', files[f]
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
                        data_all[model, var, stns[ns], dt_c] = dat_interp
                        dts_all.append(dt_c)

    #--- Get sorted unique list of all dates seen.
    dts_unique = list(sorted(set(dts_all)))
    
    #--- Create pickle file.
    pickle.dump((data_all, dts_unique, models, vars, stns), \
                open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
# Get seasonal precipitation totals and make a time series plot of them.
#---------------------------------------------------------------------------
(spr_pcp_tot, sum_pcp_tot, fal_pcp_tot, win_pcp_tot, years) = \
          get_seasonal_stats(data_all, dts_unique, models, stns, 'PREC', 'tot')

mod = 'gfdl-cm3'
plotfname = plot_dir + '/' + mod + '_season_precip.png'
titlein = mod + ', ' + str(years[0]) + ' - ' + str(years[len(years)-1])
iret = ts(win_pcp_tot, spr_pcp_tot, sum_pcp_tot, fal_pcp_tot, years, \
          stns, mod, titlein, plotfname, 'PREC')

#---------------------------------------------------------------------------
# Get seasonal temperature averages.
#---------------------------------------------------------------------------
(spr_tmx_avg, sum_tmx_avg, fal_tmx_avg, win_tmx_avg, years) = \
          get_seasonal_stats(data_all, dts_unique, models, stns, 'T2MAX', 'avg')
(spr_tmn_avg, sum_tmn_avg, fal_tmn_avg, win_tmn_avg, years) = \
          get_seasonal_stats(data_all, dts_unique, models, stns, 'T2MIN', 'avg')

(spr_spdavg_avg, sum_spdavg_avg, fal_spdavg_avg, win_spdavg_avg, years) = \
                 get_seasonal_stats(data_all, dts_unique, models, stns, \
                                    'SPDUV10MEAN', 'avg')

mod = 'gfdl-cm3'
plotfname = plot_dir + '/' + mod + '_season_tmx_avg.png'
titlein = mod + ', TMAX Average, ' + str(years[0]) + ' - ' + \
          str(years[len(years)-1])
iret = ts(win_tmx_avg, spr_tmx_avg, sum_tmx_avg, fal_tmx_avg, years, \
          stns, mod, titlein, plotfname, 'T2MAX')

mod = 'gfdl-cm3'
plotfname = plot_dir + '/' + mod + '_season_tmn_avg.png'
titlein = mod + ', TMIN Average, ' + str(years[0]) + ' - ' + \
          str(years[len(years)-1])
iret = ts(win_tmn_avg, spr_tmn_avg, sum_tmn_avg, fal_tmn_avg, years, \
          stns, mod, titlein, plotfname, 'T2MIN')

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------

sys.exit()

