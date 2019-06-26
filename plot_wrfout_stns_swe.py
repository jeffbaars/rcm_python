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
import csv

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
data_dir   = rcm_dir + '/data'
snotel_dir = data_dir + '/snotel'
run_dirs   = rcm_dir + '/rundirs'
geo_em     = rcm_dir + '/data' + '/geo_em.d02.nc'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']
station_file = data_dir + '/station_file_swe.txt'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
varswe = 'SNOW'
stat = 'tot'

sdt = '197001'
edt = '209912'

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']
mod_cols = ['indigo', 'blue', 'deepskyblue', \
            'darkgreen', 'lime', 'yellow', \
            'magenta', 'red', 'salmon', 'gray', 'darkgray', 'lightblue']

yyyy_s = sdt[0:4]
yyyy_e = edt[0:4]        

#---------------------------------------------------------------------------
# Load stations file.
#---------------------------------------------------------------------------
print 'Loading ', station_file
stns      = []
latpts    = []
lonpts    = []
elevs     = []
stn_names = []
with open(station_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        stns.append(row[0])
        latpts.append(float(row[1]))
        lonpts.append(float(row[2]))
        elevs.append(float(row[3]))

#---------------------------------------------------------------------------
# Load geo_em file.
#---------------------------------------------------------------------------
print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, read in Snow Water Equivalent data for all stns.
#---------------------------------------------------------------------------
pfswe = pickle_dir + '/swe_stns.pkl'
if os.path.isfile(pfswe):
    print 'Loading ', pfswe
    (stns, dts_all, swe_all) = pickle.load(open(pfswe, 'rb'))
else:
    for m in range(len(models)):
        model = models[m]

        #--- Get list of extracted files for this model.
        files = []
        for d in range(len(data_dirs)):
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

        #--- Initialize our big SWE array.  
        if 'swe_all' not in locals():
            swe_all = np.ones((len(models), len(stns), len(files))) * np.nan

        #--- Loop over files found above reading in Apr 1 SWE for every station.
        dts_c = []    
        for f in range(len(files)):
            print files[f]
            (vardat, ntimes) = load_wrfout_swe(files[f], [varswe])
            dt_c = re.findall(r'\.(\d{10})\.', files[f])
            dts_c.append(dt_c)

            for s in range(len(stns)):
                dat_interp = bilinear_interpolate(lat, lon, \
                                                  vardat[varswe][:,:,0], \
                                                  latpts[s], lonpts[s], 0)
                print stns[s], dat_interp
                swe_all[m,s,f] = dat_interp

        #--- Initialize our dates array if this is the first model.  Otherwise
        #--- do a date check to make sure # of dates are identical to first
        #--- model's.  Bomb otherwise, as that means a model is missing a file.
        if 'dts_all' not in locals():
            dts_all = dts_c
        else:
            if sorted(dts_c) != sorted(dts_all):
                print 'dts_c   = ', dts_c
                print 'dts_all = ', dts_all                
                sys.exit('dts problem!')

    pickle.dump((stns, dts_all, swe_all), open(pfswe,'wb'), -1)

#---------------------------------------------------------------------------
# nan out 1970, because that's a partial winter given models started 01/01/70.
#---------------------------------------------------------------------------
swe_all[:,:,0] = np.nan

#---------------------------------------------------------------------------
# Get years from dts_all.
#---------------------------------------------------------------------------
years = []
for d in range(len(dts_all)):
    years.append(dts_all[d][0][0:4])

#---------------------------------------------------------------------------
# Make plots using ts_swe.
#---------------------------------------------------------------------------
(nm,ns,nf) = swe_all.shape
for s in range(len(stns)):
    stn = stns[s]
    plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_swe.png'
    titlein = station_name_dict[stn] + ' (' + stn.upper() + '), ' + \
              labels[varswe+stat] + ', ' + \
              str(years[0]) + ' - ' + str(years[len(years)-1])

    #--- need to get obs and read them...
    obs = []
    obs_years = []
    
    iret = ts_swe(swe_all, years, obs, obs_years, stn, s, models, titlein, \
                  plotfname, varswe, stat, mod_cols)
        

sys.exit()

