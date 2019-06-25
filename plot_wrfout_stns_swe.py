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
vars = ['SNOW']

sdts = ['197001', '203001', '207001']
edts = ['200012', '206012', '209912']

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']
#models = ['ccsm4']

#---------------------------------------------------------------------------
# Make run directory, load geo_em file.
#---------------------------------------------------------------------------
rundir = run_dirs + '/plot_wrfout_grid_avg_snowrad.' + str(os.getpid())
print 'making ', rundir
os.makedirs(rundir)

#rundir = '/home/disk/spock/jbaars/rcm/rundirs/plot_wrfout_grid_avg.4149'

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, make a maps short-wave down radiation (SWDOWN).
#---------------------------------------------------------------------------
#var = 'SWDOWN'
#for m in range(len(models)):
#    model = models[m]
#    for n in range(len(sdts)):
#        sdt = sdts[n]
#        edt = edts[n]
#        yyyy_s = sdt[0:4]
#        yyyy_e = edt[0:4]        
#
#        #--- Get list of extracted files for this model.
#        files = get_extract_files(data_dirs, model, sdt, edt, 1)
#
#        #--- Get files broken out by year and season.
#        (spr_files, sum_files, fal_files, win_files, years) = \
#                    get_seasonal_files(files)
#
#        for season in seasons:
#            file_plot = rundir + '/' + model + '_' + var + '_' + season + \
#                        '_' + yyyy_s + '_' + yyyy_e + '.nc'
#            syscom = 'ncra'
#            nyears = 0
#            for y in range(len(years)):
#                yyyy = years[y]
#
#                if season == 'summer':
#                    if (sum_files[yyyy] == sum_files[yyyy]):
#                        nyears += 1
#                    for f in range(len(sum_files[yyyy])):
#                        syscom = syscom + ' ' + sum_files[yyyy][f]
#                    else:
#                        continue
#                elif season == 'spring':
#                    if (spr_files[yyyy] == spr_files[yyyy]):
#                        nyears += 1
#                        for f in range(len(spr_files[yyyy])):
#                            syscom = syscom + ' ' + spr_files[yyyy][f]
#                        else:
#                            continue
#                elif season == 'winter':
#                    if (win_files[yyyy] == win_files[yyyy]):
#                        nyears += 1
#                        for f in range(len(win_files[yyyy])):
#                            syscom = syscom + ' ' + win_files[yyyy][f]
#                        else:
#                            continue
#                elif season == 'fall':                        
#                    if (fal_files[yyyy] == fal_files[yyyy]):
#                        nyears += 1
#                        for f in range(len(fal_files[yyyy])):
#                            syscom = syscom + ' ' + fal_files[yyyy][f]
#                        else:
#                            continue
#
#            syscom = syscom + ' ' + file_plot
#            print 'syscom = ', syscom        
#            os.system(syscom)
#
#            print 'Loading ', var, ' from ', file_plot
#            nc = NetCDFFile(file_plot, 'r')
#            ntimes = len(nc.dimensions['Time'])
#            vardat = {}
#            dat_tmp = nc.variables[var][:,:,:]
#            dat_tmp = np.transpose(dat_tmp)
#            vardat[var] = dat_tmp
#            nc.close()
#            
#            titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
#                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
#    
#            plotdat = vardat[var][:,:,0]
#            if (var == 'PREC'):
#                plotdat = plotdat * mm2in
#    
#            for z in range(len(zooms)):
#                plotfname = rundir + '/' + model + '_' + sdt + '_' + \
#                            edt + '_' + season + '_' + var + '_' + \
#                            zooms[z] + '.png'
#                mapper(var, lat, lon, plotdat, levs_swdown, cmap_swdown, \
#                       maplw, colorbar_labs[var], titlein, plotfname, \
#                       zooms[z])
##            sys.exit()

'''
#---------------------------------------------------------------------------
# For each model, make a maps Snow Water Equivalent.
#---------------------------------------------------------------------------
var = 'SNOW'
for m in range(len(models)):

    model = models[m]

    for n in range(len(sdts)):
        sdt = sdts[n]
        edt = edts[n]
        yyyy_s = sdt[0:4]
        yyyy_e = edt[0:4]        

        #--- Get list of extracted files for this model.
        files = []
        for d in range(len(data_dirs)):
#            path_c = data_dirs[d] + '/' + model + '/extract/*0401*.nc'
            path_c = data_dirs[d] + '/' + model + '/swe/*0401*.nc'
            files_c = glob.glob(path_c)
        
            for f in range(len(files_c)):
                if 'snowrad' not in files_c[f]:
                    continue

                dt_c = re.findall(r'\.(\d{10})\.', files_c[f])
                if (dt_c):
                    dt_c = dt_c[0]
                    yyyymm_c = dt_c[0:6]
                    if (yyyymm_c >= sdt and yyyymm_c <= edt):
                        files.append(files_c[f])

        #--- Make plot for each zooms.
        file_plot = rundir + '/' + model + '_' + var + '_apr1_' + \
                    yyyy_s + '_' + yyyy_e + '.nc'
        syscom = 'ncra'
        for f in range(len(files)):
            syscom = syscom + ' ' + files[f]
            
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
                  'Apr 1st, ' + yyyy_s + ' - ' + yyyy_e
    
        plotdat = vardat[var][:,:,0]
    
        for z in range(len(zooms)):
            plotfname = rundir + '/' + model + '_' + sdt + '_' + \
                        edt + '_apr1_' + var + '_' + \
                        zooms[z] + '.png'
            mapper(var, lat, lon, plotdat, levs_swe, cmap_swe, \
                   maplw, colorbar_labs[var], titlein, plotfname, \
                   zooms[z])
'''
                  
#---------------------------------------------------------------------------
# Make map of Snow Water Equivalent for ensemble mean.
#---------------------------------------------------------------------------
#(nx,ny) = lat.shape
#snowrad_mean = np.ones((nx, ny, len(sdts))) * np.nan
#var = 'SNOW'
#for n in range(len(sdts)):
#    sdt = sdts[n]
#    edt = edts[n]
#    yyyy_s = sdt[0:4]
#    yyyy_e = edt[0:4]        
#
#    files = []
#    dts_all = []
#    for m in range(len(models)):
#        model = models[m]
#        for d in range(len(data_dirs)):
#            path_c = data_dirs[d] + '/' + model + '/swe/*0401*.nc'
#            files_c = glob.glob(path_c)
#        
#            for f in range(len(files_c)):
#                if 'snowrad' not in files_c[f]:
#                    continue
#
#                dt_c = re.findall(r'\.(\d{10})\.', files_c[f])
#
#                if (dt_c):
#                    dt_c = dt_c[0]
#                    yyyymm_c = dt_c[0:6]
#                    if (yyyymm_c >= sdt and yyyymm_c <= edt):
#                        files.append(files_c[f])
#                        dts_all.append(dt_c)
#                        
#    #--- Make plot for each zooms.
##    file_plot = rundir + '/ensmean_' + var + '_apr1_' + \
##                yyyy_s + '_' + yyyy_e + '.nc'
#
#    for f in range(len(files)):
#        print files[f], dts_all[f]
#        (vardat, ntimes) = load_wrfout_swe(files[f], [var])
#        plotdat = vardat[var][:,:,0]
#        nx,ny = plotdat.shape
#        if f == 0:
#            dat_all = np.ones((nx, ny, len(files))) * np.nan
#        dat_all[:,:,f] = plotdat
#        print 'np.amax(plotdat) = ', np.amax(plotdat)
#
#    snowrad_mean[:,:,n] = np.mean(dat_all, axis=2)
#    snowrad_std = np.std(dat_all, axis=2)    
#
#    print 'np.amax(snowrad_std)  = ', np.amax(snowrad_std)
#    print 'np.amax(snowrad_mean) = ', np.amax(snowrad_mean[:,:,n])
#    print 'np.amax(dat_all[:,:,0]) = ', np.amax(dat_all[:,:,0])
#    print 'np.amax(dat_all[:,:,1]) = ', np.amax(dat_all[:,:,1])            
#
#    for z in range(len(zooms)):
#        titlein = 'Ensemble Mean, ' + var_lab[var] + ', ' + \
#                  'Apr 1st, ' + yyyy_s + ' - ' + yyyy_e
#        plotfname = rundir + '/ensmean_' + sdt + '_' + \
#                    edt + '_apr1_' + var + '_' + zooms[z] + '.png'
#        mapper(var, lat, lon, snowrad_mean[:,:,n], levs_swe, cmap_swe, \
#               maplw, colorbar_labs[var], titlein, plotfname, \
#               zooms[z])
#
#        titlein = 'Ensemble Standard Deviation, ' + var_lab[var] + ', ' + \
#                  'Apr 1st, ' + yyyy_s + ' - ' + yyyy_e
#        plotfname = rundir + '/ensstd_' + sdt + '_' + \
#                    edt + '_apr1_' + var + '_' + zooms[z] + '.png'
#        mapper(var, lat, lon, snowrad_std, levs_swe, cmap_swe, \
#               maplw, colorbar_labs[var], titlein, plotfname, \
#               zooms[z])
#

pf = pickle_dir + '/swe_mean.pkl'
#pickle.dump((snowrad_mean), open(pf,'wb'), -1)

(snowrad_mean) = pickle.load(open(pf, 'rb'))
var = 'SNOW'
for n in range(1,len(sdts)):
    diff_c = snowrad_mean[:,:,n] - snowrad_mean[:,:,0]
    dttit = '(' + sdts[n][0:4] + '-' + edts[n][0:4] + ') - ' + \
            '(' + sdts[0][0:4] + '-' + edts[0][0:4] + ')'

    print diff_c
    print np.amin(diff_c)
    print np.amax(diff_c)
    print diff_c.shape
    
    for z in range(len(zooms)):
        titlein = 'Apr 1st ' + var_lab[var] + \
                  ', Ensemble Mean Difference, ' + ', ' + dttit
        plotfname = rundir + '/ensmean_diff_' + sdts[n] + '_' + \
                    edts[n] + '_apr1_' + var + '_' + zooms[z] + '.png'
        mapper(var, lat, lon, diff_c, levs_swe_diff, cmap_swe_diff, \
               maplw, colorbar_labs[var], titlein, plotfname, \
               zooms[z])

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

