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
obs_dir    = rcm_dir + '/obs'
run_dirs   = rcm_dir + '/rundirs'
geo_em     = rcm_dir + '/data' + '/geo_em.d02.nc'
data_dirs  = ['/home/disk/r2d2/steed/cmip5/rcp8.5', \
              '/home/disk/vader/steed/cmip5/rcp8.5', \
              '/home/disk/jabba/steed/cmip5/rcp8.5']
station_file = obs_dir + '/station_file_swe.txt'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
varswe = 'SNOW'

sdt = '197001'
edt = '209912'

models = ['access1.0', 'access1.3', 'bcc-csm1.1', 'canesm2', \
          'ccsm4', 'csiro-mk3.6.0', 'fgoals-g2', 'gfdl-cm3', \
          'giss-e2-h', 'miroc5', 'mri-cgcm3', 'noresm1-m']

yyyy_s = sdt[0:4]
yyyy_e = edt[0:4]        

#---------------------------------------------------------------------------
# Load stations file.
#---------------------------------------------------------------------------
print 'Loading ', station_file
stns   = []
latpts = []
lonpts = []
elevs  = []
with open(station_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        stns.append(row[0])
        latpts.append(float(row[1]))
        lonpts.append(float(row[2]))
        elevs.append(float(row[3]))

##---------------------------------------------------------------------------
## Make run directory, load geo_em file.
##---------------------------------------------------------------------------
#rundir = run_dirs + '/plot_wrfout_grid_avg_snowrad.' + str(os.getpid())
#print 'making ', rundir
#os.makedirs(rundir)

print 'Loading ', geo_em
(lat, lon, hgt) = load_geo_em(geo_em)

#---------------------------------------------------------------------------
# For each model, read in Snow Water Equivalent data for all stns.
#---------------------------------------------------------------------------
pfswe = pickle_dir + '/swe_stns.pkl'
if os.path.isfile(pfswe):
    print 'Loading ', pfswe
    (stns, dts, swe_all) = pickle.load(open(pfswe, 'rb'))
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
# 
#---------------------------------------------------------------------------


sys.exit()

