#!/usr/bin/python
import  sys, os, os.path, time, glob, re, math
import numpy as np
#import matplotlib.pyplot as plt
#from make_cmap import *
from netCDF4 import Dataset as NetCDFFile
from utils_wrfout import *
from utils_load_data import *
from utils_plot import *
from utils_date import *
import pickle
from collections import defaultdict

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
#vars = ['T2MAX', 'T2MIN', 'PREC']
#vars = ['T2MAX', 'T2MIN']
vars = ['PREC']
#vars = ['SPDUV10MAX', 'SPDUV10MEAN']

sdts = ['197001', '203001', '207001']
edts = ['200012', '206012', '209912']

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']

dtfmt = "%Y%m%d%H"

#---------------------------------------------------------------------------
# Make run directory, load geo_em file.
#---------------------------------------------------------------------------
rundir = run_dirs + '/plot_wrfout_grid_avg.' + str(os.getpid())
print 'making ', rundir
os.makedirs(rundir)

#rundir = run_dirs + '/plot_wrfout_grid_avg.23227' #
#rundir = run_dirs + '/plot_wrfout_grid_avg.16978' # PREC
#rundir = run_dirs + '/plot_wrfout_grid_avg.27465'

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

'''
#---------------------------------------------------------------------------
# For each model, make a map of every season...
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
        (spr_files, sum_files, fal_files, win_files, ann_files, years) = \
                    get_seasonal_files(files)

        for var in vars:
            nyears_ann = 0
            for season in seasons:
                file_plot = rundir + '/' + model + '_' + var + '_' + season + \
                            '_' + yyyy_s + '_' + yyyy_e + '.nc'

                if (var == 'T2MAX' or var == 'T2MIN'):
                    syscom = 'ncra'
                elif (var == 'PREC'):
                    syscom = 'ncra -y ttl'

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
                    elif season == 'annual':                        
                        if (ann_files[yyyy] == ann_files[yyyy]):
                            nyears += 1
                            for f in range(len(ann_files[yyyy])):
                                syscom = syscom + ' ' + ann_files[yyyy][f]
                            else:
                                continue

                nyears_ann = nyears_ann + nyears

                syscom = syscom + ' ' + file_plot
                print 'syscom = ', syscom        
                os.system(syscom)
        
                print 'Loading ', var, ' from ', file_plot
                (vardat, ntimes, dt) = load_wrfout(file_plot, [var])
        
                titlein = model.upper() + ', Average ' + var_lab[var] + ', ' + \
                          seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
        
                plotdat = vardat[var][:,:,0]
                if (var == 'PREC'):
                    if 'prec_tot' in locals():
                        prec_tot = np.add(prec_tot, plotdat)
                    else:
                        prec_tot = plotdat
                    
                    plotdat = (plotdat / nyears) * mm2in
                    levs_c = levs_pcp
                    cmap_c = cmap_pcp
                else:
                    plotdat = plotdat - 273.15
                    levs_c = levs_temp
                    cmap_c = cmap_temp
        
                for z in range(len(zooms)):
                    plotfname = rundir + '/' + model + '_' + sdt + '_' + \
                                edt + '_' + season + '_' + var + '_' + \
                                zooms[z] + '.png'
                    mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
                           maplw, colorbar_labs[var], titlein, plotfname, \
                           zooms[z])
#                sys.exit()

            if (var == 'PREC'):
                file_plot = rundir + '/' + model + '_' + var + '_annual_' + \
                            yyyy_s + '_' + yyyy_e + '.nc'
                f = NetCDFFile(file_plot, 'w')
                (nx,ny) = prec_tot.shape
                f.createDimension('nx', nx)
                f.createDimension('ny', ny)
                pcpnc = f.createVariable('pcp', 'f', ('nx','ny'))
                pcpnc[:,:] = prec_tot
                pcpnc.description = 'total annual precip'
                pcpnc.units       = 'mm'
                f.close()

sys.exit()
'''

#---------------------------------------------------------------------------
# Ensemble spread plots for every season...
#---------------------------------------------------------------------------
for n in range(len(sdts)):
    spr_files = {}
    sum_files = {}
    fal_files = {}
    win_files = {}
    ann_files = {}    
    
    sdt = sdts[n]
    edt = edts[n]
    yyyy_s = sdt[0:4]
    yyyy_e = edt[0:4]        

    #--- Fill seasonal file dictionaries for each model.
    for m in range(len(models)):
        model = models[m]

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dirs, model, sdt, edt, 0)

        #--- Get files broken out by year and season.
        (spr_files[model], sum_files[model], fal_files[model], \
         win_files[model], ann_files[model], years) = get_seasonal_files(files)

    #--- For each variable and season, make a plot.
    for var in vars:
        if (var == 'T2MAX' or var == 'T2MIN'):
            levs_c = levs_temp_std
            cmap_c = cmap_temp_std
        elif (var == 'PREC'):
            levs_c = levs_pcp_std
            cmap_c = cmap_pcp_std

        for season in seasons:
            titlein = 'Ensemble Standard Deviation, ' + var_lab[var] + ', ' + \
                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e
            pickle_file = pickle_dir + '/map_all_std_' + yyyy_s + '_' + \
                          yyyy_e + '_' + season + '_' + var + '.pickle'
            if os.path.isfile(pickle_file):
                print 'loading ', pickle_file
                (varmean) = pickle.load(open(pickle_file, 'rb'))
            else:
                (varmean) = get_ensmean_std(var, models, years, season, \
                                            sum_files, spr_files, \
                                            win_files, fal_files, \
                                            ann_files, pickle_file)
            for z in range(len(zooms)):
                plotfname = rundir + '/map_all_std_' + sdt + '_' + edt + '_' + \
                            season + '_' + var + '_' + zooms[z] + '.png'
                mapper(var, lat, lon, varmean, levs_c, cmap_c, \
                       maplw, colorbar_labs[var], titlein, plotfname, \
                       zooms[z])
sys.exit()

'''
#---------------------------------------------------------------------------
# Ensemble mean plots for every season...
#---------------------------------------------------------------------------
for n in range(len(sdts)):
    spr_files = {}
    sum_files = {}
    fal_files = {}
    win_files = {}
    ann_files = {}    
    
    sdt = sdts[n]
    edt = edts[n]
    yyyy_s = sdt[0:4]
    yyyy_e = edt[0:4]        

    #--- Fill seasonal file dictionaries for each model.
    for m in range(len(models)):
        model = models[m]

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dirs, model, sdt, edt, 0)

        #--- Get files broken out by year and season.
        (spr_files[model], sum_files[model], fal_files[model], \
         win_files[model], ann_files[model], years) = get_seasonal_files(files)

    #--- For each variable and season, make a plot.
    for var in vars:
        if (var == 'T2MAX' or var == 'T2MIN' or var == 'SPDUV10MAX' or \
            var == 'SPDUV10MEAN'):
            ncocom = 'ncra'
        elif (var == 'PREC'):
            ncocom = 'ncra -y ttl'

        for season in seasons:
            file_plot = rundir + '/' + var + '_' + season + \
                        '_' + yyyy_s + '_' + yyyy_e + '.nc'

            (syscom, nyears) = get_ensmean_syscom(file_plot, var, models, \
                                                  years, season, \
                                                  sum_files, spr_files, \
                                                  win_files, fal_files, \
                                                  ann_files, ncocom)
            if os.path.isfile(file_plot):
                print 'file_plot already exists.  file_plot = ', file_plot
            else:
                iret = run_ensmean_syscom(file_plot, var, ncocom, syscom, \
                                          yyyy_s, yyyy_e, season, rundir)

            #--- Add number of years (nmodels * years) as dimension to nc file.
            syscom = 'ncap2 -s \'defdim(\"nmodyears\",' + str(nyears) + \
                ')\'' + ' ' +  file_plot
            print syscom
            os.system(syscom)
            
            print 'Loading ', var, ' from ', file_plot
            (vardat, ntimes, dt) = load_wrfout(file_plot, [var])

            titlein = 'Ensemble Mean, ' + var_lab[var] + ', ' + \
                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e

            plotdat = vardat[var][:,:,0]
            if (var == 'PREC'):
                plotdat = (plotdat / nyears) * mm2in
                levs_c = levs_pcp
                cmap_c = cmap_pcp
            elif var == 'T2MAX' or var == 'T2MIN':
                plotdat = plotdat - 273.15
                levs_c = levs_temp
                cmap_c = cmap_temp
            elif var == 'SPDUV10MAX':
                levs_c = levs_spd
                cmap_c = cmap_spd

            for z in range(len(zooms)):
                plotfname = rundir + '/map_all_' + sdt + '_' + edt + '_' + \
                            season + '_' + var + '_' + zooms[z] + '.png'

                print var, lat, lon, plotdat, levs_c, cmap_c, \
                       maplw, colorbar_labs[var], titlein, plotfname, \
                       zooms[z]
                
                mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
                       maplw, colorbar_labs[var], titlein, plotfname, \
                       zooms[z])
            #sys.exit()
#sys.exit()
'''
                
#---------------------------------------------------------------------------
# Get difference plots for the two sdts / edts periods.
#---------------------------------------------------------------------------
for nd in range(1,len(sdts)):
    print nd, sdts[nd]
    sdt_hist = sdts[0]
    edt_hist = edts[0]
    yyyy_s_hist = sdt_hist[0:4]
    yyyy_e_hist = edt_hist[0:4]        
    sdt_fore = sdts[nd]
    edt_fore = edts[nd]
    yyyy_s_fore = sdt_fore[0:4]
    yyyy_e_fore = edt_fore[0:4]        
    
    for var in vars:
        for season in seasons:

            #--- Get diff between forecast and historical for var and season.
            file_hist = rundir + '/' + var + '_' + season + \
                        '_' + yyyy_s_hist + '_' + yyyy_e_hist + '.nc'
            forename = yyyy_s_fore + '_' + yyyy_e_fore
            file_fore = rundir + '/' + var + '_' + season + \
                        '_' + forename + '.nc'
            file_diff = rundir + '/' + var + '_' + season + '_' + \
                        forename + '_diff.nc'
            syscom = 'ncdiff -O ' + file_fore + ' ' + file_hist + ' ' +file_diff
            print 'syscom = ', syscom
            os.system(syscom)

            #--- Load diff file.
            print 'Loading ', var, ' from ', file_diff
            (vardat, ntimes, dt) = load_wrfout(file_diff, [var])

            plotdat = vardat[var][:,:,0]

            nameadd = ''
            if (var == 'PREC'):
                (vardat, ntimes, dt) = load_wrfout(file_hist, [var])
                vardat_hist = vardat[var][:,:,0]
                nmodyears = get_nc_dim(file_hist, 'nmodyears')
                vardat_hist = (vardat_hist / nmodyears) * mm2in

                (vardat, ntimes, dt) = load_wrfout(file_fore, [var])
                vardat_fore = vardat[var][:,:,0]
                nmodyears = get_nc_dim(file_fore, 'nmodyears')
                vardat_fore = (vardat_fore / nmodyears) * mm2in

                #--- Do plot of percentage difference:
                plotdat = ((vardat_fore - vardat_hist) / \
                           ((vardat_fore + vardat_hist) / 2)) * 100
                levs_c = levs_pcp_percdiff
                cmap_c = cmap_pcp_percdiff
                var_long = 'Total Precipitation'
                cb_lab = 'Total Precip Percentage Difference (%)'
                stat_c = 'Percentage Difference'
                units_c = '%'
                nameadd = '_PERCDIFF'

#                #--- Do plot of actual difference:
#                plotdat = vardat_fore - vardat_hist
#                cb_lab = colorbar_labs[var]
#                if season == 'annual' or season == 'fall' or season == 'winter':
#                    levs_c = levs_pcp_diff_annual
#                    cmap_c = cmap_pcp_diff_annual
#                else:
#                    levs_c = levs_pcp_diff
#                    cmap_c = cmap_pcp_diff
#                var_long = 'Total Precipitation'
#                stat_c = 'Difference'
#                units_c = 'in.'
                
            elif var == 'T2MAX':
                cb_lab = colorbar_labs[var]
                levs_c = levs_temp_diff
                cmap_c = cmap_temp_diff
                var_long = 'Average Maximum Temperature'
                stat_c = 'Difference'
                units_c = '$^\circ$C'
            elif var == 'T2MIN':            
                cb_lab = colorbar_labs[var]
                levs_c = levs_temp_diff
                cmap_c = cmap_temp_diff
                var_long = 'Average Minimum Temperature'
                stat_c = 'Difference'
                units_c = '$^\circ$C'
            elif var == 'SPDUV10MAX':            
                diff_c = plotdat
                (vardat, ntimes, dt) = load_wrfout(file_hist, [var])
                vardat_hist = vardat[var][:,:,0]
                (vardat, ntimes, dt) = load_wrfout(file_fore, [var])
                vardat_fore = vardat[var][:,:,0]            
                plotdat = ((vardat_fore - vardat_hist) / \
                           ((vardat_fore + vardat_hist) / 2)) * 100
                levs_c = levs_spd_diff
                cmap_c = cmap_spd_diff
                var_long = 'Average Maximum Wind Speed'
                cb_lab = colorbar_labs[var]
                stat_c = 'Difference'
                units_c = 'm/s'

#            titlein = 'Ensemble Mean, ' + var_lab[var] + ', ' + \
#                      seasons_lab[season] + ', ' + yyyy_s + ' - ' + yyyy_e

            titlein =  seasons_lab[season] + ' ' + var_long + ' ' + stat_c + \
                      ' (' + units_c + ')' + \
                      ', (' + yyyy_s_fore + ' - ' + yyyy_e_fore + ') ' + \
                      ' - ' + \
                      '(' + yyyy_s_hist + ' - ' + yyyy_e_hist + ')'
            
            for z in range(len(zooms)):
                plotfname = rundir + '/map_diff_' + sdt_fore + '_' + \
                            edt_fore + '_' + season + '_' + var + nameadd + \
                            '_' + zooms[z] + '.png'
                mapper(var, lat, lon, plotdat, levs_c, cmap_c, \
                       maplw, cb_lab, titlein, plotfname, zooms[z])
            #sys.exit()
sys.exit()

