#!/usr/bin/python
import  sys, os, os.path, time, glob, re, math
import numpy as np
#import matplotlib.pyplot as plt
from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_date import *
import pickle
from collections import defaultdict
from make_cmap import *

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
run_dirs   = rcm_dir + '/rundirs'
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

#---------------------------------------------------------------------------
# Make a color map using make_cmap.
#---------------------------------------------------------------------------
#if (varname == 'tos'):
#    levs_c = levs_temp
#    cmap_c = cmap_temp
#elif(varname == 'ua'):
#    levs_c = levs_prec
#    cmap_c = cmap_prec
levs_c = levs_temp
cmap_c = cmap_temp
        
levs_norm_c = []
for i in range(len(levs_c)):
    x = float(levs_c[i])
    norm_c = (x - min(levs_c)) / (max(levs_c) - min(levs_c))
    levs_norm_c.append(norm_c)

my_cmap = make_cmap(cmap_c, bit = True, position = levs_norm_c)

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
    data_dir_c = data_dir + '/' + model + '/extract'

    #--- Get list of extracted files for this model.
    files = get_extract_files(data_dir_c, model, sdt, edt)

    #--- Get files broken out by year and season.
    (spr_files, sum_files, fal_files, win_files, years) = \
                get_seasonal_files(files)

    #--- Loop over year, making maps of each 
    for y in range(len(years)):
        yyyy = years[y]

        #--- Winter plots.
        if (win_files[yyyy] == win_files[yyyy]):
            test = run_mapper('T2MAX', model, rundir, 'winter', win_files, \
                              yyyy, 'Temperature (K)', lat, lon, levs_c, \
                              my_cmap)
            
        if (sum_files[yyyy] == sum_files[yyyy]):
            test = run_mapper('T2MAX', model, rundir, 'summer', sum_files, \
                              yyyy, 'Temperature (K)', lat, lon, levs_c, \
                              my_cmap)

            test = run_mapper('PREC', model, rundir, 'summer', sum_files, \
                              yyyy, 'Precipitation (in.)', lat, lon, levs_c, \
                              my_cmap)
            sys.exit()

        sys.exit()
            
#            file_plot = rundir + '/sum_' + yyyy + '.nc'
#            syscom = 'ncra'
#            for f in range(len(sum_files[yyyy])):
#                print yyyy, sum_files[yyyy][f]
#                syscom = syscom + ' ' + sum_files[yyyy][f]
#            syscom = syscom + ' ' + file_plot
#            print 'syscom = ', syscom
#            os.system(syscom)
#
#            var = 'T2MAX'
#            print 'Loading ', var, ' from ', file_plot
#            (vardat, ntimes, dt) = load_wrfout(file_plot, [var])
#
#            colorbar_lab = 'Temperature (K)'
#            titlein = model.upper() + \
#                      ', Average Max Temperature (K), Summer (JJA) ' + yyyy
#
#            plotdat = vardat[var][:,:,0]
#            for z in range(len(zooms)):
#                plotfname = rundir + '/sum_' + var + '_' + yyyy + '_' + \
#                            zooms[z] + '.png'
#                mapper(lat, lon, plotdat, levs_c, my_cmap, maplw, \
#                       colorbar_lab, titlein, plotfname, zooms[z])

    sys.exit()


sys.exit()

