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
(lat, lon) = load_geo_em(geo_em)

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
            file_plot = rundir + '/win_' + yyyy + '.nc'
            syscom = 'ncra'
            for f in range(len(win_files[yyyy])):
                print yyyy, win_files[yyyy][f]
                syscom = syscom + ' ' + win_files[yyyy][f]
            syscom = syscom + ' ' + file_plot
            print 'syscom = ', syscom
            os.system(syscom)

            var = 'T2MAX'
            print 'Loading ', var, ' from ', file_plot
            (vardat, ntimes, dt) = load_wrfout(file_plot, [var])

            colorbar_lab = 'Temperature (K)'
            titlein = model + ', Average Max Temperature (K), DJF ' + yyyy
            plotfname = rundir + '/win_' + var + '_' + yyyy + '.png'

            plotdat = vardat[var][:,:,0]
            mapper(lat, lon, plotdat, minlat, maxlat, minlon, maxlon, \
                   levs_c, my_cmap, maplw, colorbar_lab, titlein, plotfname)

        if (sum_files[yyyy] == sum_files[yyyy]):
            file_plot = rundir + '/sum_' + yyyy + '.nc'
            syscom = 'ncra'
            for f in range(len(sum_files[yyyy])):
                print yyyy, sum_files[yyyy][f]
                syscom = syscom + ' ' + sum_files[yyyy][f]
            syscom = syscom + ' ' + file_plot
            print 'syscom = ', syscom
            os.system(syscom)

            var = 'T2MAX'
            print 'Loading ', var, ' from ', file_plot
            (vardat, ntimes, dt) = load_wrfout(file_plot, [var])

            colorbar_lab = 'Temperature (K)'
            titlein = model + ', Average Max Temperature (K), JJA ' + yyyy
            plotfname = rundir + '/sum_' + var + '_' + yyyy + '.png'

            plotdat = vardat[var][:,:,0]
            mapper(lat, lon, plotdat, minlat, maxlat, minlon, maxlon, \
                   levs_c, my_cmap, maplw, colorbar_lab, titlein, plotfname)



    sys.exit()



#    #--- Loop over all extract files, reading them in.
#    springs = defaultdict(list)
#    summers = defaultdict(list)
#    falls   = defaultdict(list)
#    winters = defaultdict(list)
#    spring_ndays = {}
#    summer_ndays = {}
#    fall_ndays   = {}
#    winter_ndays = {}
#
#    for f in range(len(files)):
##        print files[f]
#        dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]
#        if (dt_c):
#            yyyy = dt_c[0:4]
#            mm   = dt_c[4:6]
#            years_mod_stn.append(yyyy)
#
#            print dt_c, yyyy, mm
#            if (mm == '03' or mm == '04' or mm == '05'):
#                if ((yyyy) in spring_ndays):
#                    spring_ndays[yyyy] += 1
#                else:
#                    spring_ndays[yyyy] = 1
#                springs[yyyy].append(files[f])
#            elif (mm == '06' or mm == '07' or mm == '08'):
#                if ((yyyy) in summer_ndays):
#                    summer_ndays[yyyy] += 1
#                else:
#                    summer_ndays[yyyy] = 1
#                summers[yyyy].append(files[f])
#            elif (mm == '09' or mm == '10' or mm == '11'):
#                if ((yyyy) in fall_ndays):
#                    fall_ndays[yyyy] += 1
#                else:
#                    fall_ndays[yyyy] = 1
#                falls[yyyy].append(files[f])
#            elif (mm == '12'):
#                if ((yyyyw) in winter_ndays):
#                    winter_ndays[yyyyw] += 1
#                else:
#                    winter_ndays[yyyyw] = 1
#                yyyyw = str(int(yyyy) + 1)
#                winters[yyyyw].append(files[f])
#            elif (mm == '01' or mm == '02'):
#                if ((yyyy) in winter_ndays):
#                    winter_ndays[yyyy] += 1
#                else:
#                    winter_ndays[yyyy] = 1
#                winters[yyyy].append(files[f])
#                
#    print springs
#    print 'springs[2020] = ', springs['2020']
#    print 'summers[2020] = ', summers['2020']
#    print 'winters[2020] = ', winters['2020']
#
#    sys.exit()

##---------------------------------------------------------------------------
## Loop over each models, grabbing data from extract wrfout files.
##---------------------------------------------------------------------------
#pickle_file = pickle_dir + '/extract_' + sdt + '_' + edt + '.pickle'
#if os.path.isfile(pickle_file):
#    (data_all, dts_unique, models, vars, stns) = \
#               pickle.load(open(pickle_file, 'rb'))
#else:
#    print 'Loading ', geo_em
#    (lat, lon) = load_geo_em(geo_em)
#
#    data_all = {}
#    dts_all = []
#    for m in range(len(models)):
#
#        model = models[m]
#        data_dir_c = data_dir + '/' + model + '/extract'
#
#        #--- Get list of extracted files for this model.
#        files = get_extract_files(data_dir_c, model, sdt, edt)
#
#        #--- Loop over all extract files, reading them in.
#        for f in range(len(files)):
#
#            print 'Loading ', files[f]
#            (vardat, ntimes, dt) = load_wrfout(files[f], vars)
#
#            #--- Loop over each vars, interpolating data to stns for
#            #--- each ntimes.
#            for v in range(len(vars)):
#
#                var = vars[v]
#                vardat_c = vardat[var]
#                
#                #--- Loop over each ntime, grabbing interpolated data for
#                #--- variable var for each stns.
#                for n in range(ntimes):
#                    dt_c = time_increment(dt, n, dtfmt)
#                    for ns in range(len(stns)):
#                        dat_interp = bilinear_interpolate(lat, lon, \
#                                                          vardat_c[:,:,n],\
#                                                          latpts[ns], \
#                                                          lonpts[ns], 0)
#                        data_all[model, var, stns[ns], dt_c] = dat_interp
#                        dts_all.append(dt_c)
#
#    #--- Get sorted unique list of all dates seen.
#    dts_unique = list(sorted(set(dts_all)))
#    
#    #--- Create pickle file.
#    pickle.dump((data_all, dts_unique, models, vars, stns), \
#                open(pickle_file,'wb'), -1)
#
##---------------------------------------------------------------------------
## Get seasonal precipitation totals and make a time series plot of them.
##---------------------------------------------------------------------------
#(spr_pcp_tot, sum_pcp_tot, fal_pcp_tot, win_pcp_tot, years) = \
#          get_seasonal_stats(data_all, dts_unique, models, stns, 'PREC', 'tot')
#
#mod = 'gfdl-cm3'
#plotfname = plot_dir + '/' + mod + '_season_precip.png'
#titlein = mod + ', ' + str(years[0]) + ' - ' + str(years[len(years)-1])
#iret = ts(win_pcp_tot, spr_pcp_tot, sum_pcp_tot, fal_pcp_tot, years, \
#          stns, mod, titlein, plotfname, 'PREC')
#
##---------------------------------------------------------------------------
## Get seasonal temperature averages.
##---------------------------------------------------------------------------
#(spr_tmx_avg, sum_tmx_avg, fal_tmx_avg, win_tmx_avg, years) = \
#          get_seasonal_stats(data_all, dts_unique, models, stns, 'T2MAX', 'avg')
#(spr_tmn_avg, sum_tmn_avg, fal_tmn_avg, win_tmn_avg, years) = \
#          get_seasonal_stats(data_all, dts_unique, models, stns, 'T2MIN', 'avg')
#
#(spr_spdavg_avg, sum_spdavg_avg, fal_spdavg_avg, win_spdavg_avg, years) = \
#                 get_seasonal_stats(data_all, dts_unique, models, stns, \
#                                    'SPDUV10MEAN', 'avg')
#
#mod = 'gfdl-cm3'
#plotfname = plot_dir + '/' + mod + '_season_tmx_avg.png'
#titlein = mod + ', TMAX Average, ' + str(years[0]) + ' - ' + \
#          str(years[len(years)-1])
#iret = ts(win_tmx_avg, spr_tmx_avg, sum_tmx_avg, fal_tmx_avg, years, \
#          stns, mod, titlein, plotfname, 'T2MAX')
#
#mod = 'gfdl-cm3'
#plotfname = plot_dir + '/' + mod + '_season_tmn_avg.png'
#titlein = mod + ', TMIN Average, ' + str(years[0]) + ' - ' + \
#          str(years[len(years)-1])
#iret = ts(win_tmn_avg, spr_tmn_avg, sum_tmn_avg, fal_tmn_avg, years, \
#          stns, mod, titlein, plotfname, 'T2MIN')

sys.exit()

