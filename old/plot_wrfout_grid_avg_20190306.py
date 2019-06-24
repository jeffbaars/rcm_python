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

##print cmap_pcpd
#print cmap_pcp
#
#sys.exit()

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
#vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#vars = ['PREC', 'T2MAX', 'T2MIN']
vars = ['T2MAX', 'T2MIN']
#vars = ['PREC']
#vars = ['SWDOWN']
#vars = ['PREC']

sdts = ['197001', '207001']
edts = ['200012', '209912']
#sdts = ['207001']
#edts = ['209912']

models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', \
            'access1.3', 'canesm2', 'noresm1-m', \
            'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2']

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
for i in range(len(levs_temp_diff)):
    x = float(levs_temp_diff[i])
    norm_c = (x - min(levs_temp_diff)) / \
             (max(levs_temp_diff) - min(levs_temp_diff))
    levs_norm_c.append(norm_c)
cmap_temp_diff = make_cmap(cmap_temp_diff, bit = True, position = levs_norm_c)

levs_norm_c = []
for i in range(len(levs_pcp)):
    x = float(levs_pcp[i])
    norm_c = (x - min(levs_pcp)) / (max(levs_pcp) - min(levs_pcp))
    levs_norm_c.append(norm_c)
cmap_pcp = make_cmap(cmap_pcp, bit = True, position = levs_norm_c)

#---------------------------------------------------------------------------
# Make run directory, load geo_em file.
#---------------------------------------------------------------------------
#rundir = run_dirs + '/plot_wrfout_grid_avg.' + str(os.getpid())
#print 'making ', rundir
#os.makedirs(rundir)

#rundir = run_dirs + '/plot_wrfout_grid_avg.17103'
#rundir = run_dirs + '/plot_wrfout_grid_avg.19749'
#rundir = run_dirs + '/plot_wrfout_grid_avg.23178'
rundir = run_dirs + '/plot_wrfout_grid_avg.16504'

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

##---------------------------------------------------------------------------
## For each model, make a maps of every season...
##---------------------------------------------------------------------------
#for m in range(len(models)):
#    model = models[m]
#    for n in range(len(sdts)):
#        sdt = sdts[n]
#        edt = edts[n]
#        yyyy_s = sdt[0:4]
#        yyyy_e = edt[0:4]        
#
#        #--- Get list of extracted files for this model.
#        files = get_extract_files(data_dirs, model, sdt, edt, 0)
#
#        #--- Get files broken out by year and season.
#        (spr_files, sum_files, fal_files, win_files, years) = \
#                    get_seasonal_files(files)
#       
#        for var in vars:
#            for season in seasons:
#                file_plot = rundir + '/' + model + '_' + var + '_' + season + \
#                            '_' + yyyy_s + '_' + yyyy_e + '.nc'
#
#                if (var == 'T2MAX' or var == 'T2MIN'):
#                    syscom = 'ncra'
#                elif (var == 'PREC'):
#                    syscom = 'ncra -y ttl'
#
#                nyears = 0
#                for y in range(len(years)):
#                    yyyy = years[y]
#                    if season == 'summer':
#                        if (sum_files[yyyy] == sum_files[yyyy]):
#                            nyears += 1
#                            for f in range(len(sum_files[yyyy])):
#                                syscom = syscom + ' ' + sum_files[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'spring':
#                        if (spr_files[yyyy] == spr_files[yyyy]):
#                            nyears += 1
#                            for f in range(len(spr_files[yyyy])):
#                                syscom = syscom + ' ' + spr_files[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'winter':
#                        if (win_files[yyyy] == win_files[yyyy]):
#                            nyears += 1
#                            for f in range(len(win_files[yyyy])):
#                                syscom = syscom + ' ' + win_files[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'fall':                        
#                        if (fal_files[yyyy] == fal_files[yyyy]):
#                            nyears += 1
#                            for f in range(len(fal_files[yyyy])):
#                                syscom = syscom + ' ' + fal_files[yyyy][f]
#                            else:
#                                continue
#                        
#                syscom = syscom + ' ' + file_plot
#                print 'syscom = ', syscom        
#                os.system(syscom)
#        
#                print 'Loading ', var, ' from ', file_plot
#                (vardat, ntimes, dt) = load_wrfout(file_plot, [var])
#        
#                titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
#                          seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
#        
#                plotdat = vardat[var][:,:,0]
#                if (var == 'PREC'):
#                    plotdat = (plotdat / nyears) * mm2in
#                    levs_c = levs_pcp
#                    cmap_c = cmap_pcp
#                else:
#                    levs_c = levs_temp
#                    cmap_c = cmap_temp
#        
#                for z in range(len(zooms)):
#                    plotfname = rundir + '/' + model + '_' + sdt + '_' + \
#                                edt + '_' + season + '_' + var + '_' + \
#                                zooms[z] + '.png'
#                    mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
#                           maplw, colorbar_labs[var], titlein, plotfname, \
#                           zooms[z])

#---------------------------------------------------------------------------
# Average all models for every season...
#---------------------------------------------------------------------------
#for n in range(len(sdts)):
#    spr_files = {}
#    sum_files = {}
#    fal_files = {}
#    win_files = {}
#    
#    sdt = sdts[n]
#    edt = edts[n]
#    yyyy_s = sdt[0:4]
#    yyyy_e = edt[0:4]        
#
#    for m in range(len(models)):
#        model = models[m]
#
#        #--- Get list of extracted files for this model.
#        files = get_extract_files(data_dirs, model, sdt, edt, 0)
#
#        #--- Get files broken out by year and season.
#        (spr_files[model], sum_files[model], fal_files[model], \
#         win_files[model], years) = get_seasonal_files(files)
#
#    for var in vars:
#        for season in seasons:
#
#            file_plot = rundir + '/' + var + '_' + season + \
#                        '_' + yyyy_s + '_' + yyyy_e + '.nc'
#
#            if (var == 'T2MAX' or var == 'T2MIN'):
#                syscom = 'ncra'
#            elif (var == 'PREC'):
#                syscom = 'ncra -y ttl'
#
#            for m in range(len(models)):
#                model = models[m]
#                sum_files_c = sum_files[model]
#                spr_files_c = spr_files[model]
#                win_files_c = win_files[model]
#                fal_files_c = fal_files[model]                
#                
#                nyears = 0
#                for y in range(len(years)):
#                    yyyy = years[y]
#                    if season == 'summer':
#                        if (sum_files_c[yyyy] == sum_files_c[yyyy]):
#                            nyears += 1
#                            for f in range(len(sum_files_c[yyyy])):
#                                syscom = syscom + ' ' + sum_files_c[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'spring':
#                        if (spr_files_c[yyyy] == spr_files_c[yyyy]):
#                            nyears += 1
#                            for f in range(len(spr_files_c[yyyy])):
#                                syscom = syscom + ' ' + spr_files_c[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'winter':
#                        if (win_files_c[yyyy] == win_files_c[yyyy]):
#                            nyears += 1
#                            for f in range(len(win_files_c[yyyy])):
#                                syscom = syscom + ' ' + win_files_c[yyyy][f]
#                            else:
#                                continue
#                    elif season == 'fall':                        
#                        if (fal_files_c[yyyy] == fal_files_c[yyyy]):
#                            nyears += 1
#                            for f in range(len(fal_files_c[yyyy])):
#                                syscom = syscom + ' ' + fal_files_c[yyyy][f]
#                            else:
#                                continue
#                        
#            syscom = syscom + ' ' + file_plot
#            print 'syscom = ', syscom
#            os.system(syscom)
#
##            #--- FOR TESTING ONLY!  REMOVE!
##            nyears = 30
#        
#            print 'Loading ', var, ' from ', file_plot
#            (vardat, ntimes, dt) = load_wrfout(file_plot, [var])
#        
#            titlein = var_lab[var] + ', ' + \
#                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
#        
#            plotdat = vardat[var][:,:,0]
#            if (var == 'PREC'):
#                plotdat = (plotdat / nyears) * mm2in
#                levs_c = levs_pcp
#                cmap_c = cmap_pcp
#            else:
#                plotdat = plotdat - 273.15
#                levs_c = levs_temp
#                cmap_c = cmap_temp
#    
#            for z in range(len(zooms)):
#                plotfname = rundir + '/map_all_' + sdt + '_' + edt + '_' + \
#                            season + '_' + var + '_' + zooms[z] + '.png'
#                mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
#                       maplw, colorbar_labs[var], titlein, plotfname, \
#                       zooms[z])
##                sys.exit()
#
#sys.exit()

#---------------------------------------------------------------------------
# Get difference plots for the two sdts / edts periods.
#---------------------------------------------------------------------------
sdt_hist = sdts[0]
edt_hist = edts[0]
yyyy_s_hist = sdt_hist[0:4]
yyyy_e_hist = edt_hist[0:4]        
sdt_fore = sdts[1]
edt_fore = edts[1]
yyyy_s_fore = sdt_fore[0:4]
yyyy_e_fore = edt_fore[0:4]        

for var in vars:
    for season in seasons:

        #--- Get diff between forecast and historical for this var and season.
        file_hist = rundir + '/' + var + '_' + season + \
                    '_' + yyyy_s_hist + '_' + yyyy_e_hist + '.nc'
        file_fore = rundir + '/' + var + '_' + season + \
                    '_' + yyyy_s_fore + '_' + yyyy_e_fore + '.nc'
        file_diff = rundir + '/' + var + '_' + season + '_' + 'diff.nc'

        syscom = 'ncdiff ' + file_fore + ' ' + file_hist + ' ' + file_diff
        print 'syscom = ', syscom
        os.system(syscom)

        #--- Load diff file.
        print 'Loading ', var, ' from ', file_diff
        (vardat, ntimes, dt) = load_wrfout(file_diff, [var])

        plotdat = vardat[var][:,:,0]        
        
        if (var == 'PREC'):
            diff_c = plotdat
            (vardat, ntimes, dt) = load_wrfout(file_hist, [var])
            vardat_hist = vardat[var][:,:,0]
            (vardat, ntimes, dt) = load_wrfout(file_fore, [var])
            vardat_fore = vardat[var][:,:,0]            
            plotdat = ((vardat_fore - vardat_hist) / \
                   ((vardat_fore + vardat_hist) / 2)) * 100

#            plotdat = plotdat * mm2in
#            levs_c = levs_pcp_diff
#            cmap_c = cmap_pcp_diff
            levs_c = levs_pcp
            cmap_c = cmap_pcp
            cb_lab = 'Total Precip Percentage Difference (%)'
            title_supp = 'Percentage Difference'
            title_supp2 = 'vs.'            
        else:
            cb_lab = colorbar_labs[var]
            levs_c = levs_temp_diff
            cmap_c = cmap_temp_diff
            title_supp = 'Difference'
            title_supp2 = '-'

#        titlein = var_lab[var] + ' ' + title_supp + ', ' + \
        titlein =  seasons_lab[season] + ' Total Precipitation ' + \
                  title_supp + ' (%)' + \
                  ', (' + yyyy_s_fore + ' - ' + yyyy_e_fore + ') ' + \
                  + ' vs. ' + \
                  '(' + yyyy_s_hist + ' - ' + yyyy_e_hist + ')'

        for z in range(len(zooms)):
            plotfname = rundir + '/map_diff_' + \
                        season + '_' + var + '_' + zooms[z] + '.png'
            mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
                   maplw, cb_lab, titlein, plotfname, zooms[z])
#        sys.exit()

sys.exit()

