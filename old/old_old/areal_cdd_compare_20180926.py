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
stat = 'avg'  #-- only option for CDD data...
sdt = '197001'
edt = '201712'
#edt = '198012'

#vars_mod = ['PREC', 'T2MAX', 'T2MIN']
vars_mod = ['T2MAX', 'T2MIN']
#vars_mod = ['PREC']
models   = ['gfdl-cm3', 'miroc5', 'bcc-csm1.1', 'access1.3', 'canesm2',
            'noresm1-m', 'ccsm4']
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
cdd_data = {}
for var in vars_mod:
    for s in range(len(seasons)):
        sm = seasons_mo[s]
        season = seasons[s]
        cdd_file = data_dir + '/cdd_' + var + '_' + sm + '_WA.dat'
        print 'reading ', cdd_file
        cdd = open(cdd_file)
        csvReader = csv.reader(cdd)
        yyyy = []
        val_avg = []
        for row in csvReader:
            (yyyy_c,val) = row[0].split()
            yyyy.append(int(yyyy_c))
            if (float(val) == -9.8999996 or float(val) == -999.98999):
                val_avg.append(np.nan)
            else:
                val_avg.append(float(val))        
        cdd.close()
        cdd_yyyy[var,season] = yyyy
        cdd_data[var,season] = val_avg

#---------------------------------------------------------------------------
# For each model, read in grids and keep "slice" of it as defined by
# mask_file read in above.
#---------------------------------------------------------------------------
sums_all = {}
cnts_all = {}
years_all = []
for m in range(len(models)):
    mod = models[m]
    files = get_extract_files(data_dirs, mod, sdt, edt)            
    files = list(sorted(set(files)))
    (seafiles, years) = get_seasonal_files_new(files)

    years_all.append(years)

    for y in range(len(years)):
        yyyy = years[y]
        print 'loading ', mod, yyyy
        for s in range(len(seasons)):
            season = seasons[s]
            if (seafiles[season,yyyy] == seafiles[season,yyyy]):
                (sums, cnts) = sum_wrfout_slice(seafiles[season,yyyy], \
                                                vars_mod, maski, maskj)
                for var in vars_mod:
                    key = (mod,season,var,yyyy)

                    print 'key = ', key
                    print 'sums, cnts = ', sums, cnts
                    sys.exit()
                    
                    #--- Winter of first year will only have Jan-Feb.  So
                    #--- nan it out.
                    if (y == 0 and season == 'winter'):
                        sums_all[key] = np.nan
                        cnts_all[key] = np.nan
                    else:
                        if (var in sums):
                            sums_all[key] = sums[var]
                            cnts_all[key] = cnts[var]

years_all = list(sorted(set(years_all[0])))
pickle.dump((sums_all,cnts_all,years_all), open(pickle_file,'wb'), -1)

sys.exit()

if os.path.isfile(pickle_file):
    print 'reading pickle file ', pickle_file
    (sums_all,cnts_all,years_all) = pickle.load(open(pickle_file, 'rb'))

#key = (mod,season,var,str(yyyy))
#
#print sums_all[key]
#sys.exit()
#for var in vars_mod:
#    if (var == 'PREC'):
#        print 'hi'
#sys.exit()


#season = 'spring'
#key = ('gfdl-cm3', season, 'PREC', '1984')
#print season, sums_all[key], cnts_all[key], sums_all[key]/cnts_all[key]
#season = 'fall'
#key = ('gfdl-cm3', season, 'PREC', '1984')
#print season, sums_all[key], cnts_all[key], sums_all[key]/cnts_all[key]
#
#sys.exit()

for s in range(len(seasons)):
    season = seasons[s]
    key = ('gfdl-cm3', season, 'PREC', '1986')
    print season, sums_all[key], cnts_all[key], sums_all[key]/cnts_all[key]
    print season, sums_all[key] * mm2in, cnts_all[key], \
          (sums_all[key]/cnts_all[key]) * mm2in

#for key, value in sums_all.items() :
#    print (key, value)

sys.exit()

#---------------------------------------------------------------------------
# Make plot for each seasons.
#---------------------------------------------------------------------------
lw_sm = 1.8
smooth_fact = 5     #-- odd number, used for smoothing time series plots.
for var in vars_mod:
    for s in range(len(seasons)):

        fig, ax = plt.subplots( figsize=(10,6) )

        season = seasons[s]
        cdd_yyyy_c = cdd_yyyy[var,season]
    
        if (var == 'T2MAX' or var == 'T2MIN'):
            cdd_c = (np.asarray(cdd_data[var,season]) - 32.0) * (5./9.)
        elif (var == 'PREC'):
            cdd_c = np.asarray(cdd_data[var,season])
    
#        plt.plot(cdd_c, alpha=0.4, label='Observed', \
#                 color='black', linestyle='None', marker='o', markersize=3)
        lab_c = 'Observed ' + str(smooth_fact) + '-pt smoothing'
        plt.plot(smooth(cdd_c, smooth_fact), label=lab_c, color='black', \
                 linewidth=lw_sm)

        for m in range(len(models)):
            mod = models[m]
            modp = []
            xlabs = []
            for y in range(len(cdd_yyyy_c)):
                yyyy = cdd_yyyy_c[y]
                xlabs.append(yyyy)
                key = (mod,season,var,str(yyyy))
                if key in sums_all:
                    if (var == 'T2MAX' or var == 'T2MIN'):                    
                        tc = (sums_all[key] / cnts_all[key]) - 273.15
                        modp.append(tc)
                    elif (var == 'PREC'):
                        tc = (sums_all[key] / cnts_all[key]) * mm2in
                        modp.append(tc)                        
                else:
                    modp.append(np.nan)

#            plt.plot(modp, alpha = 0.17, label = mod.upper(), \
#                     color = mod_cols[m], linewidth=lw_sm)
            lab_c = mod.upper() + ', ' + str(smooth_fact) + '-pt smoothing'
            plt.plot(smooth(modp, smooth_fact), label = lab_c, \
                     color = mod_cols[m], linewidth=lw_sm)

        #--- y-axis labeling.
        plt.ylim(ylims[var+stat])
        plt.ylabel('Seasonal ' + labels[var+stat], fontsize=fs+1)
        plt.tick_params(axis='y', which='major', labelsize=fs+1)    

        #--- x-axis labels.
        xticks_c = range(0,len(xlabs), 1)
        xlabs_c = [ xlabs[i] for i in xticks_c ]    
        plt.xticks(xticks_c)
        plt.tick_params(axis='x', which='major', labelsize=fs-1)
        ax.set_xticklabels(xlabs_c, rotation=90)        
        plt.xlabel('Year', fontsize=fs+1)

        titlein = 'Climate Division Data Areal Verification, WA state, ' + \
                  var + ', ' + seasons_lab[season] + ', ' + \
                  labels[var+stat] + ', ' + sdt + ' - ' + edt
        plt.title(titlein, fontsize=titlefs-1, fontweight='bold')

        plt.tight_layout()
        plt.grid()
        plt.legend(fontsize = fs-1, loc = 'best')

        plotfname = plot_dir + '/areal_' + var + '_' + season + '.png'
        print 'xli ', plotfname, ' &'
        plt.savefig(plotfname)
        plt.close()

#        sys.exit()
sys.exit()

