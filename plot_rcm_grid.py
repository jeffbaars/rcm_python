#!/usr/bin/python
import  sys, os, os.path, time, glob, re, math
import numpy as np
#import matplotlib.pyplot as plt
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
from utils_cmap import *
import pickle
from collections import defaultdict
from make_cmap import *
import matplotlib as mpl
mpl.use('Agg')

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
run_dirs   = rcm_dir + '/rundirs'
geo_em     = rcm_dir + '/data' + '/geo_em.d02.nc'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
vars = ['PREC', 'T2MAX', 'T2MIN']

sdt = '197001'
edt = '209912'

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']

#---------------------------------------------------------------------------
# Make run directory, load geo_em file.
#---------------------------------------------------------------------------
rundir = run_dirs + '/plot_wrfout_grid.' + str(os.getpid())
print 'making ', rundir
os.makedirs(rundir)

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, make a maps of every season...
#---------------------------------------------------------------------------
for m in range(len(models)):
    model = models[m]

    #--- Get list of extracted files for this model.
    files = get_extract_files(data_dirs, model, sdt, edt, 0)

    #--- Get files broken out by year and season.
    (spr_files, sum_files, fal_files, win_files, ann_files, years) = \
                get_seasonal_files(files)

    #--- Loop over year, making maps of each 
    for y in range(len(years)):
        yyyy = years[y]

        #--- Winter plots.
        if (win_files[yyyy] == win_files[yyyy]):
            iret = run_mapper('T2MAX', model, rundir, 'winter', win_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('T2MIN', model, rundir, 'winter', win_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('PREC', model, rundir, 'winter', win_files,\
                              yyyy, lat, lon, levs_pcp, cmap_pcp)

        #--- Summer plots.
        if (sum_files[yyyy] == sum_files[yyyy]):
            iret = run_mapper('T2MAX', model, rundir, 'summer', sum_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('T2MIN', model, rundir, 'summer', sum_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('PREC', model, rundir, 'summer', sum_files,\
                              yyyy, lat, lon, levs_pcp, cmap_pcp)

        #--- Spring plots.
        if (spr_files[yyyy] == spr_files[yyyy]):
            iret = run_mapper('T2MAX', model, rundir, 'spring', spr_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('T2MIN', model, rundir, 'spring', spr_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('PREC', model, rundir, 'spring', spr_files,\
                              yyyy, lat, lon, levs_pcp, cmap_pcp)

        #--- Fall plots.
        if (fal_files[yyyy] == fal_files[yyyy]):
            iret = run_mapper('T2MAX', model, rundir, 'fall', fal_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('T2MIN', model, rundir, 'fall', fal_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('PREC', model, rundir, 'fall', fal_files,\
                              yyyy, lat, lon, levs_pcp, cmap_pcp)

        #--- Annual plots.
        if (ann_files[yyyy] == ann_files[yyyy]):
            iret = run_mapper('T2MAX', model, rundir, 'annual', ann_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('T2MIN', model, rundir, 'annual', ann_files,\
                              yyyy, lat, lon, levs_temp, cmap_temp)
            iret = run_mapper('PREC', model, rundir, 'annual', ann_files,\
                              yyyy, lat, lon, levs_pcp, cmap_pcp)


sys.exit()

