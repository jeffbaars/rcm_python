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
geo_em     = rcm_dir + '/data' + '/geo_em.d02.nc'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#vars = ['PREC', 'T2MAX', 'T2MIN']
#vars = ['T2MAX', 'T2MIN']
vars = ['PREC']

sdts = ['197001', '207001']
edts = ['200012', '209912']
#sdts = ['207001']
#edts = ['209912']

#models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', 'access1.3']
#models = ['gfdl-cm3']
models = ['noresm1-m']
dtfmt = "%Y%m%d%H"

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
rundir = run_dirs + '/plot_wrfout_grid_avg.' + str(os.getpid())
print 'making ', rundir
os.makedirs(rundir)

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, make a maps of every season...
#---------------------------------------------------------------------------
for m in range(len(models)):
    model = models[m]
    for n in range(len(sdts)):
        sdt = sdts[n]
        edt = edts[n]
        yyyy_s = sdt[0:4]
        yyyy_e = edt[0:4]        

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dirs, model, sdt, edt, 0)

        #--- Get files broken out by year and season.
        (spr_files, sum_files, fal_files, win_files, years) = \
                    get_seasonal_files(files)
       
        for var in vars:
            for season in seasons:
                file_plot = rundir + '/' + model + '_' + var + '_' + season + \
                            '_' + yyyy_s + '_' + yyyy_e + '.nc'
                syscom = 'ncra'
                for y in range(len(years)):
                    yyyy = years[y]
                    if (sum_files[yyyy] == sum_files[yyyy]):
                        for f in range(len(sum_files[yyyy])):
                            syscom = syscom + ' ' + sum_files[yyyy][f]
                    else:
                        continue
    
                syscom = syscom + ' ' + file_plot
                print 'syscom = ', syscom        
                os.system(syscom)
        
                print 'Loading ', var, ' from ', file_plot
                (vardat, ntimes, dt) = load_wrfout(file_plot, [var])
        
                titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
                          seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
        
                plotdat = vardat[var][:,:,0]

                print plotdat
                sys.exit()
                
                if (var == 'PREC'):
                    plotdat = plotdat * mm2in
                    levs_in = levs_pcp
                    cmap_in = cmap_pcp
                else:
                    levs_in = levs_temp
                    cmap_in = cmap_temp
        
                for z in range(len(zooms)):
                    plotfname = rundir + '/' + model + '_' + sdt + '_' + \
                                edt + '_' + season + '_' + var + '_' + \
                                zooms[z] + '.png'
                    mapper(var, lat, lon, plotdat, levs_in, cmap_in, \
                           maplw, colorbar_labs[var], titlein, plotfname, \
                           zooms[z])

                    sys.exit()
       
#    if (var == 'T2MAX' or var == 'T2MIN'):
#        if (os.path.exists(file_plot)):
#            print 'nc file ', file_plot, ' already exists-- using it'
#        else:
#            syscom = 'ncra'
#            for f in range(len(files[yyyy])):
#                print yyyy, files[yyyy][f]
#                syscom = syscom + ' ' + files[yyyy][f]
#            syscom = syscom + ' ' + file_plot
#            print 'syscom = ', syscom        
#            os.system(syscom)
#    elif (var == 'PREC'):
#        if (os.path.exists(file_plot)):
#            print 'nc file ', file_plot, ' already exists-- using it'
#        else:
#            syscom = 'ncra -y ttl'
#            for f in range(len(files[yyyy])):
#                print yyyy, files[yyyy][f]
#                syscom = syscom + ' ' + files[yyyy][f]
#            syscom = syscom + ' ' + file_plot
#            print 'syscom = ', syscom        
#            os.system(syscom)


#    #--- Loop over year, making maps of each 
#    for y in range(len(years)):
#        yyyy = years[y]
#
#        #--- Winter plots.
#        if (win_files[yyyy] == win_files[yyyy]):
#            iret = run_mapper('T2MAX', model, rundir, 'winter', win_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('T2MIN', model, rundir, 'winter', win_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('PREC', model, rundir, 'winter', win_files,\
#                              yyyy, lat, lon, levs_pcp, cmap_pcp)
#
#        #--- Summer plots.
#        if (sum_files[yyyy] == sum_files[yyyy]):
#            iret = run_mapper('T2MAX', model, rundir, 'summer', sum_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('T2MIN', model, rundir, 'summer', sum_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('PREC', model, rundir, 'summer', sum_files,\
#                              yyyy, lat, lon, levs_pcp, cmap_pcp)
#
#        #--- Spring plots.
#        if (spr_files[yyyy] == spr_files[yyyy]):
#            iret = run_mapper('T2MAX', model, rundir, 'spring', spr_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('T2MIN', model, rundir, 'spring', spr_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('PREC', model, rundir, 'spring', spr_files,\
#                              yyyy, lat, lon, levs_pcp, cmap_pcp)
#
#        #--- Fall plots.
#        if (fal_files[yyyy] == fal_files[yyyy]):
#            iret = run_mapper('T2MAX', model, rundir, 'fall', fal_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('T2MIN', model, rundir, 'fall', fal_files,\
#                              yyyy, lat, lon, levs_temp, cmap_temp)
#            iret = run_mapper('PREC', model, rundir, 'fall', fal_files,\
#                              yyyy, lat, lon, levs_pcp, cmap_pcp)



sys.exit()

