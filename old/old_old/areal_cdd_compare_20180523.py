#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import pickle
from utils_date import *
from utils_wrfout import *
import csv

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
data_dir   = rcm_dir + '/data'
cdd_data   = data_dir + '/mask_pts.dat'
mask_file  = data_dir + '/mask_pts.dat'

geo_em    = '/home/disk/a125/steed/run/geo_em.d02.nc'

data_dirs = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']

pickle_file = pickle_dir + '/T2MAX_slice.pickle'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
sdt = '197001'
#edt = '201712'
edt = '197512'

#vars_mod = ['PREC', 'T2MAX', 'T2MIN']
vars_mod = ['T2MAX']
models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', 'access1.3']
mod_cols = ['b',        'r',      'g',          'y']

#---------------------------------------------------------------------------
# Read mask lat/lon points file.
#---------------------------------------------------------------------------
mask = open(mask_file)
csvReader = csv.reader(mask)
maski = []
maskj = []
for row in csvReader:
    maski.append(int(row[0]))
    maskj.append(int(row[1]))
mask.close()

#---------------------------------------------------------------------------
# Read Climate Division Data file.
#---------------------------------------------------------------------------
cdd_yyyy = {}
cdd_tmx_avg = {}
for s in range(len(seasons)):
    sm = seasons_mo[s]
    season = seasons[s]
    cdd_file = data_dir + '/cdd_' + sm + '_WA.dat'
    cdd = open(cdd_file)
    csvReader = csv.reader(cdd)
    yyyy = []
    tmx_avg = []
    for row in csvReader:
        (yyyy_c,val) = row[0].split()
        yyyy.append(int(yyyy_c))
        tmx_avg.append(float(val))
    cdd.close()
    cdd_yyyy[season] = yyyy
    cdd_tmx_avg[season] = tmx_avg
    
##---------------------------------------------------------------------------
## Load geo_em file.
##---------------------------------------------------------------------------
#print 'Loading ', geo_em
#nc = NetCDFFile(geo_em, 'r')
#lat = nc.variables['XLAT_M'][0,:,:]
#lon = nc.variables['XLONG_M'][0,:,:]
#hgt = nc.variables['HGT_M'][0,:,:]    
#lat = np.transpose(lat)
#lon = np.transpose(lon)
#hgt = np.transpose(hgt)
#nc.close()

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
sums_all = {}
cnts_all = {}
years_all = []
for m in range(len(models)):
    mod = models[m]

    print 'Loading slice for model ', mod

    files = get_extract_files(data_dirs, mod, sdt, edt)            
    files = list(sorted(set(files)))

    (seafiles, years) = get_seasonal_files_new(files)

    years_all.append(years)

    for y in range(len(years)):
        yyyy = years[y]
        print 'loading year ', yyyy
        for s in range(len(seasons)):
            season = seasons[s]
            if (seafiles[season,yyyy] == seafiles[season,yyyy]):
                (sums, cnts) = sum_wrfout_slice(seafiles[season,yyyy], \
                                                vars_mod, maski, maskj)
                for var in vars_mod:
                    key = (mod,season,var,yyyy)
                    if (var in sums):
                        sums_all[key] = sums[var]
                        cnts_all[key] = cnts[var]

years_all = list(sorted(set(years_all[0])))
pickle.dump((sums_all,cnts_all,years_all), open(pickle_file,'wb'), -1)

sys.exit()

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
var = 'T2MAX'
stat = 'avg'
fig, ax = plt.subplots( figsize=(10,6) )
lw_sm = 1.8
smooth_fact     = 5       #-- odd number, used for smoothing time series plots.

for s in range(len(seasons)):
    season = seasons[s]
    
    cdd_yyyy_c = cdd_yyyy[season]
    cdd_tmx_avg_c = (cdd_tmx_avg[season] - 32.0) * (5./9.)

    plt.plot(cdd_tmx_avg_c, alpha=0.6, label='Observed', color='black', \
             linestyle='None', marker='o', markersize=3)
    lab_c = 'Observed ' + str(smooth_fact) + '-pt smoothing'
    plt.plot(smooth(cdd_tmx_avg_c, smooth_fact), label=lab_c, color='black', \
             linewidth=lw_sm)

    for m in range(len(models)):
        mod = models[m]
        modp = []
        xlabs = []
        for y in range(len(cdd_yyyy_c)):
            yyyy = cdd_yyyy_c[y]
            xlabs.append(yyyy)
            key = (mod,season,var,yyyy)
            if key in sums_all:
                tc = sums_all[key] / cnts_all[key]
                tf = (tc * 9./5.) + 32.0
                modp.append(tf)
            else:
                modp.append(np.nan)

        plt.plot(modp, alpha = 0.25, label = mod.upper(), color = mod_cols[m], \
                 linewidth=lw_sm)
        lab_c = mod.upper() + ', ' + str(smooth_fact) + '-pt smoothing'
        plt.plot(smooth(modp, smooth_fact), label = lab_c, color = mod_cols[m],\
                 linewidth=lw_sm)

    #--- y-axis labeling.
    plt.ylim(ylims[15,90])
    plt.ylabel('Seasonal ' + labels[var+stat], fontsize=fs+1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    

    #--- x-axis labels.
    xticks_c = range(0,len(xlabs), 4)
    xlabs_c = [ xlabs[i] for i in xticks_c ]    
    plt.xticks(xticks_c)
    plt.tick_params(axis='x', which='major', labelsize=fs-1)
    ax.set_xticklabels(xlabs_c, rotation=90)        
    plt.xlabel('Year', fontsize=fs+1)

    titlein = 'test'
    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    plt.tight_layout()
    plt.grid()
    plt.legend(fontsize = fs-1, loc = 'best')

    plotfname = './test.png'
    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)
    plt.close()


    sys.exit()

'''
    #--- Loop over year...
    for y in range(len(years)):
#        sums = {}
#        cnts = {}
        yyyy = years[y]

        print 'loading year ', yyyy
        
        if (win_files[yyyy] == win_files[yyyy]):
            season = 'winter'
#            (sums, cnts) = sum_wrfout_slice(win_files[yyyy], vars_mod, \
#                                            maski, maskj, sums, cnts)
            (sums, cnts) = sum_wrfout_slice(win_files[yyyy], vars_mod, \
                                            maski, maskj)
            for var in vars_mod:
                key = (mod,season,var,yyyy)
                sums_all[key] = sums[var]
                cnts_all[key] = cnts[var]
        if (spr_files[yyyy] == spr_files[yyyy]):
            season = 'spring'
#            (sums, cnts) = sum_wrfout_slice(spr_files[yyyy], vars_mod, \
#                                            maski, maskj, sums, cnts)
            (sums, cnts) = sum_wrfout_slice(spr_files[yyyy], vars_mod, \
                                            maski, maskj)
            for var in vars_mod:
                key = (mod,season,var,yyyy)
                sums_all[key] = sums[var]
                cnts_all[key] = cnts[var]
        if (sum_files[yyyy] == sum_files[yyyy]):
            season = 'summer'
#            (sums, cnts) = sum_wrfout_slice(sum_files[yyyy], vars_mod, \
#                                            maski, maskj, sums, cnts)
            (sums, cnts) = sum_wrfout_slice(sum_files[yyyy], vars_mod, \
                                            maski, maskj)
            for var in vars_mod:
                key = (mod,season,var,yyyy)
                sums_all[key] = sums[var]
                cnts_all[key] = cnts[var]
        if (fal_files[yyyy] == fal_files[yyyy]):
            season = 'fall'
#            (sums, cnts) = sum_wrfout_slice(fal_files[yyyy], vars_mod, \
#                                            maski, maskj, sums, cnts)
            (sums, cnts) = sum_wrfout_slice(fal_files[yyyy], vars_mod, \
                                            maski, maskj)
            for var in vars_mod:
                key = (mod,season,var,yyyy)
                sums_all[key] = sums[var]
                cnts_all[key] = cnts[var]

years_all = list(sorted(set(years_all[0])))
pickle.dump((sums_all,cnts_all,years_all), open(pickle_file,'wb'), -1)
sys.exit()
'''
            
#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
(sums_all,cnts_all,years_all) = pickle.load(open(pickle_file, 'rb'))

fig, ax = plt.subplots( figsize=(10,6) )

plotfname = './test.png'

season = 'summer'
var = 'T2MAX'
for m in range(len(models)):
    mod = models[m]
    vals = []
    for yyyy in years_all:
        key = (mod,season,var,yyyy)
        if key in sums_all:
            avg_c = sums_all[key] / float(cnts_all[key])
            if (var == 'PREC'):
                avg_pr = avg_c * mm2in
            elif (var == 'T2MAX' or var == 'T2MIN'):
                avg_pr = ((avg_c - 273.15) * (9./5.)) + 32.0
                print mod, yyyy, season, var, avg_pr
        else:
            avg_pr = np.nan
        vals.append(avg_pr)

    plt.plot(vals, label = mod.upper(), color = mod_cols[m], linewidth=2)

#    sys.exit()

plt.tight_layout()
plt.grid()
plt.legend(fontsize = fs-1, loc = 'best')

print 'xli ', plotfname, ' &'
plt.savefig(plotfname)
plt.close()

sys.exit()

