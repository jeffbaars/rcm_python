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
var = 'PREC'
sdt = '202001'
edt = '202912'
models = ['gfdl-cm3']
dtfmt = "%Y%m%d%H"

stns = ['KSEA', 'KSMP', 'KYKM']
latpts = [47.44472, 47.27667, 46.56417]
lonpts = [-122.31361, -121.33722, -120.53361]
#stns = ['KSEA']
#latpts = [47.44472]
#lonpts = [-122.31361]
#stns = ['KSMP']
#latpts = [47.27667]
#lonpts = [-121.33722]

#---------------------------------------------------------------------------
# Loop over each models, grabbing data from extract wrfout files.
#---------------------------------------------------------------------------
pickle_file = pickle_dir + '/' + var + '_' + sdt + '_' + edt + '.pickle'
if os.path.isfile(pickle_file):
    (data_all, dts_unique, models, stns) = pickle.load(open(pickle_file, 'rb'))
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
            (vardat, ntimes, dt) = load_wrfout(files[f], var)

            #--- Loop over each ntime, grabbing interpolated data for
            #--- variable var for each stns.
            for n in range(ntimes):
                dt_c = time_increment(dt, n, dtfmt)
                for ns in range(len(stns)):
                    prec_interp = bilinear_interpolate(lat, lon, vardat[:,:,n],\
                                                       latpts[ns], lonpts[ns],\
                                                       0)
                    data_all[model, stns[ns], dt_c] = prec_interp
                    dts_all.append(dt_c)

    #--- Get sorted unique list of all dates seen.
    dts_unique = list(sorted(set(dts_all)))
    
    #--- Create pickle file.
    pickle.dump((data_all, dts_unique, models, stns),open(pickle_file,'wb'), -1)

#---------------------------------------------------------------------------
# Get seasonal totals.
#---------------------------------------------------------------------------
(springs, summers, falls, winters, years) = \
          get_seasonal_totals(data_all, dts_unique, models, stns)

#---------------------------------------------------------------------------
# Plot times series of total precip for each station.
#---------------------------------------------------------------------------
mod = 'gfdl-cm3'
plotfname = 'test.png'
titlein = mod + ', ' + str(years[0]) + ' - ' + str(years[len(years)-1])
iret = ts(winters, springs, summers, falls, years, stns, mod, \
          titlein, plotfname)

#---------------------------------------------------------------------------
# 
#---------------------------------------------------------------------------

sys.exit()

