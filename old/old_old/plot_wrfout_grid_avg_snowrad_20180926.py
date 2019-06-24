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
#vars = ['SNOW', 'SWDOWN']
#vars = ['SNOW']
vars = ['SWDOWN']

sdts = ['197001', '207001']
edts = ['200012', '209912']
#sdts = ['207001']
#edts = ['209912']

#models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', 'access1.3']
#models   = ['gfdl-cm3', 'miroc5', 'access1.3']
#models = ['gfdl-cm3']
models = ['access1.3']
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

levs_norm_c = []
for i in range(len(levs_swe)):
    x = float(levs_swe[i])
    norm_c = (x - min(levs_swe)) / (max(levs_swe) - min(levs_swe))
    levs_norm_c.append(norm_c)
cmap_swe = make_cmap(cmap_swe, bit = True, position = levs_norm_c)

levs_norm_c = []
for i in range(len(levs_swdown)):
    x = float(levs_swdown[i])
    norm_c = (x - min(levs_swdown)) / (max(levs_swdown) - min(levs_swdown))
    levs_norm_c.append(norm_c)
cmap_swdown = make_cmap(cmap_swdown, bit = True, position = levs_norm_c)

#---------------------------------------------------------------------------
# Make run directory, load geo_em file.
#---------------------------------------------------------------------------
rundir = run_dirs + '/plot_wrfout_grid_avg.' + str(os.getpid())
print 'making ', rundir
os.makedirs(rundir)

#rundir = '/home/disk/spock/jbaars/rcm/rundirs/plot_wrfout_grid_avg.4149'

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, make a maps short-wave down radiation (SWDOWN).
#---------------------------------------------------------------------------
var = 'SWDOWN'
for m in range(len(models)):
    model = models[m]
    for n in range(len(sdts)):
        sdt = sdts[n]
        edt = edts[n]
        yyyy_s = sdt[0:4]
        yyyy_e = edt[0:4]        

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dirs, model, sdt, edt, 1)

        #--- Get files broken out by year and season.
        (spr_files, sum_files, fal_files, win_files, years) = \
                    get_seasonal_files(files)
       
        for season in seasons:
            file_plot = rundir + '/' + model + '_' + var + '_' + season + \
                        '_' + yyyy_s + '_' + yyyy_e + '.nc'
            syscom = 'ncra'
            nyears = 0
            for y in range(len(years)):
                yyyy = years[y]

                if season == 'summer':
                    if (sum_files[yyyy] == sum_files[yyyy]):
                        nyears += 1
                    for f in range(len(sum_files[yyyy])):
                        syscom = syscom + ' ' + sum_files[yyyy][f]
                    else:
                        continue
                elif season == 'spring':
                    if (spr_files[yyyy] == spr_files[yyyy]):
                        nyears += 1
                        for f in range(len(spr_files[yyyy])):
                            syscom = syscom + ' ' + spr_files[yyyy][f]
                        else:
                            continue
                elif season == 'winter':
                    if (win_files[yyyy] == win_files[yyyy]):
                        nyears += 1
                        for f in range(len(win_files[yyyy])):
                            syscom = syscom + ' ' + win_files[yyyy][f]
                        else:
                            continue
                elif season == 'fall':                        
                    if (fal_files[yyyy] == fal_files[yyyy]):
                        nyears += 1
                        for f in range(len(fal_files[yyyy])):
                            syscom = syscom + ' ' + fal_files[yyyy][f]
                        else:
                            continue

            syscom = syscom + ' ' + file_plot
            print 'syscom = ', syscom        
            os.system(syscom)

            print 'Loading ', var, ' from ', file_plot
            nc = NetCDFFile(file_plot, 'r')
            ntimes = len(nc.dimensions['Time'])
            vardat = {}
            dat_tmp = nc.variables[var][:,:,:]
            dat_tmp = np.transpose(dat_tmp)
            vardat[var] = dat_tmp
            nc.close()
            
            titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
    
            plotdat = vardat[var][:,:,0]
            if (var == 'PREC'):
                plotdat = plotdat * mm2in
    
            for z in range(len(zooms)):
                plotfname = rundir + '/' + model + '_' + sdt + '_' + \
                            edt + '_' + season + '_' + var + '_' + \
                            zooms[z] + '.png'
                mapper(var, lat, lon, plotdat, levs_swdown, cmap_swdown, \
                       maplw, colorbar_labs[var], titlein, plotfname, \
                       zooms[z])
#            sys.exit()

##---------------------------------------------------------------------------
## For each model, make a maps Snow Water Equivalent.
##---------------------------------------------------------------------------
#var = 'SNOW'
#for m in range(len(models)):
#
#    model = models[m]
#
#    for n in range(len(sdts)):
#        sdt = sdts[n]
#        edt = edts[n]
#        yyyy_s = sdt[0:4]
#        yyyy_e = edt[0:4]        
#
#        #--- Get list of extracted files for this model.
#        files = []
#        for d in range(len(data_dirs)):
#            path_c = data_dirs[d] + '/' + model + '/extract/*0401*.nc'
#            files_c = glob.glob(path_c)
#        
#            for f in range(len(files_c)):
#                if 'snowrad' not in files_c[f]:
#                    continue
#
#                dt_c = re.findall(r'\.(\d{10})\.', files_c[f])
#                if (dt_c):
#                    dt_c = dt_c[0]
#                    yyyymm_c = dt_c[0:6]
#                    if (yyyymm_c >= sdt and yyyymm_c <= edt):
#                        files.append(files_c[f])
#
#        #--- Make plot for each zooms.
#        file_plot = rundir + '/' + model + '_' + var + '_apr1_' + \
#                    yyyy_s + '_' + yyyy_e + '.nc'
#        syscom = 'ncra'
#        for f in range(len(files)):
#            syscom = syscom + ' ' + files[f]
#            
#        syscom = syscom + ' ' + file_plot
#        print 'syscom = ', syscom        
#        os.system(syscom)
#    
#        print 'Loading ', var, ' from ', file_plot
#        nc = NetCDFFile(file_plot, 'r')
#        ntimes = len(nc.dimensions['Time'])
#        vardat = {}
#        dat_tmp = nc.variables[var][:,:,:]
#        dat_tmp = np.transpose(dat_tmp)
#        vardat[var] = dat_tmp
#        nc.close()
#        titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
#                  'Apr 1st, ' + yyyy_s + ' - ' + yyyy_e
#    
#        plotdat = vardat[var][:,:,0]
#    
#        for z in range(len(zooms)):
#            plotfname = rundir + '/' + model + '_' + sdt + '_' + \
#                        edt + '_apr1_' + var + '_' + \
#                        zooms[z] + '.png'
#            mapper(var, lat, lon, plotdat, levs_swe, cmap_swe, \
#                   maplw, colorbar_labs[var], titlein, plotfname, \
#                   zooms[z])
#
##---------------------------------------------------------------------------
## Make map of difference.
##---------------------------------------------------------------------------
#for m in range(len(models)):
#    model = models[m]
#
#    syscom = 'ncdiff'
#    file_diff = rundir + '/' + model + '_' + var + '_apr1_' + 'diff.nc'
#
#    for n in range(len(sdts)):
#        sdt = sdts[n]
#        edt = edts[n]
#        yyyy_s = sdt[0:4]
#        yyyy_e = edt[0:4]        
#        file_plot = rundir + '/' + model + '_' + var + '_apr1_' + \
#                    yyyy_s + '_' + yyyy_e + '.nc'
#            if os.path.isfile(file_plot):
#                syscom = syscom + ' ' + file_plot
#
#    syscom = syscom + ' ' + file_diff
#
#    print 'syscom = ', syscom        
#    os.system(syscom)
#
#    print 'Loading ', var, ' from ', file_diff
#    nc = NetCDFFile(file_diff, 'r')
#    ntimes = len(nc.dimensions['Time'])
#    vardat = {}
#    dat_tmp = nc.variables[var][:,:,:]
#    dat_tmp = np.transpose(dat_tmp)
#    vardat[var] = dat_tmp
#    nc.close()
#
#    titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
#              'Apr 1st, Difference, ' + yyyy_s + ' - ' + yyyy_e
#        
#    plotdat = vardat[var][:,:,0]
#    
#    for z in range(len(zooms)):
#        plotfname = rundir + '/' + model + '_' + sdt + '_' + \
#                    edt + '_apr1_diff_' + var + '_' + \
#                    zooms[z] + '.png'
#        mapper(var, lat, lon, plotdat, levs_swe, cmap_swe, \
#               maplw, colorbar_labs[var], titlein, plotfname, \
#               zooms[z])

sys.exit()

