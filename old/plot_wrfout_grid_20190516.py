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
#vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
vars = ['PREC', 'T2MAX', 'T2MIN']
#vars = ['T2MIN']
#sdt = '202001'
#edt = '202912'
#sdt = '197001'
sdt = '204501'
edt = '209912'
#models = ['gfdl-cm3']
models = ['miroc5']
dtfmt = "%Y%m%d%H"

#var_c = 'T2MIN'
##print 'levs[var_c] = ', levs[var_c]
#print 'cmaps[var_c] = ', cmaps[var_c]
#test =  cmaps[var_c].extend(cmaps[var_c])
#print 'test = ', test
#sys.exit()

#---------------------------------------------------------------------------
# Make a color map using make_cmap.
#---------------------------------------------------------------------------
levs_norm_c = []
for i in range(len(levs_temp)):
    x = float(levs_temp[i])
    norm_c = (x - min(levs_temp)) / (max(levs_temp) - min(levs_temp))
    levs_norm_c.append(norm_c)
cmap_temp = make_cmap(cmap_temp, bit = True, position = levs_norm_c)

levs_norm_c = []
for i in range(len(levs_pcp)):
    x = float(levs_pcp[i])
    norm_c = (x - min(levs_pcp)) / (max(levs_pcp) - min(levs_pcp))
    levs_norm_c.append(norm_c)
cmap_pcp = make_cmap(cmap_pcp, bit = True, position = levs_norm_c)

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



#        #--- Winter plots.
#        if (win_files[yyyy] == win_files[yyyy]):
#            for v in range(len(vars)):
#                var_c = vars[v]
#                iret = run_mapper(vars[v], model, rundir, 'winter', win_files,\
#                                  yyyy, lat, lon, levs[var_c], cmaps[var_c])
#
#        #--- Summer plots.
#        if (sum_files[yyyy] == sum_files[yyyy]):
#            for v in range(len(vars)):
#                var_c = vars[v]
#                iret = run_mapper(var_c, model, rundir, 'summer', sum_files,\
#                                  yyyy, lat, lon, \
#                                  levs[var_c], cmaps[str(var_c)])
#
#        #--- Spring plots.
#        if (spr_files[yyyy] == spr_files[yyyy]):
#            for v in range(len(vars)):
#                var_c = vars[v]
#                iret = run_mapper(vars[v], model, rundir, 'spring', spr_files,\
#                                  yyyy, lat, lon, levs[var_c], cmaps[var_c])
#            
#        #--- Fall plots.
#        if (fal_files[yyyy] == fal_files[yyyy]):
#            for v in range(len(vars)):
#                var_c = vars[v]
#                iret = run_mapper(vars[v], model, rundir, 'fall', fal_files,\
#                                  yyyy, lat, lon, levs[var_c], cmaps[var_c])
#
#        sys.exit()


sys.exit()

