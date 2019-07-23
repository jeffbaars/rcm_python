#!/usr/bin/python
import os, sys, glob, re, math
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset as NetCDFFile
import pickle
from utils_date import *
from utils_stats import *
from utils_ghcnd_obs import *

#---------------------------------------------------------------------------
# Get model pickle file name
#---------------------------------------------------------------------------
def get_mod_pkl_name (pdir, model, sdt, edt):
    pf = pdir + '/extract_' + model + '_' + sdt + '_' + edt + '.pkl'
    return pf

#---------------------------------------------------------------------------
# Load geo_em file.
#---------------------------------------------------------------------------
def load_geo_em(geo_em):
    nc = NetCDFFile(geo_em, 'r')
    lat = nc.variables['XLAT_M'][0,:,:]
    lon = nc.variables['XLONG_M'][0,:,:]
    hgt = nc.variables['HGT_M'][0,:,:]    
    lat = np.transpose(lat)
    lon = np.transpose(lon)
    hgt = np.transpose(hgt)
    nc.close()
    return lat, lon, hgt

#---------------------------------------------------------------------------
# Get netcdf dimension
#---------------------------------------------------------------------------
def get_nc_dim(ncin, dimname):
    nc = NetCDFFile(ncin, 'r')
    dim = len(nc.dimensions[dimname])
    nc.close()
    return dim

#---------------------------------------------------------------------------
# Load wrfout file.
#---------------------------------------------------------------------------
def load_wrfout(wrfout, vars):
    nc = NetCDFFile(wrfout, 'r')
    ntimes = len(nc.dimensions['Time'])

    dataout = {}
    for v in range(len(vars)):
        var = vars[v]
        dat_tmp = nc.variables[var][:,:,:]
        dat_tmp = np.transpose(dat_tmp)
        dataout[var] = dat_tmp

    #--- Get date from global attributes and clean it up and trim it down.
    dt_all = nc.START_DATE
    dt = dt_all.replace('-', '')
    dt = dt.replace(':', '')
    dt = dt.replace('_', '')        
    dt = dt[0:10]
    nc.close()
    return dataout, ntimes, dt

#---------------------------------------------------------------------------
# Load SWE wrfout file.
#---------------------------------------------------------------------------
def load_wrfout_swe(wrfout, vars):
    nc = NetCDFFile(wrfout, 'r')
    ntimes = len(nc.dimensions['Time'])

    dataout = {}
    for v in range(len(vars)):
        var = vars[v]
        dat_tmp = nc.variables[var][:,:,:]
        dat_tmp = np.transpose(dat_tmp)
        dataout[var] = dat_tmp
    nc.close()
    return dataout, ntimes

#---------------------------------------------------------------------------
# Get list of extracted wrfout files between sdt_yyyymm and edt_yyyymm.
#---------------------------------------------------------------------------
def get_extract_files(data_dirs, model, sdt_yyyymm, edt_yyyymm, snowrad):
    files_out = []
    for d in range(len(data_dirs)):
        path_c = data_dirs[d] + '/' + model + '/extract/*.nc'
        files = glob.glob(path_c)
        
        for f in range(len(files)):

            if snowrad == 0:
                if 'snowrad' in files[f]:
                    continue
            else:
                if 'snowrad' not in files[f]:
                    continue

            dt_c = re.findall(r'\.(\d{10})\.', files[f])
            if (dt_c):
                dt_c = dt_c[0]
                yyyymm_c = dt_c[0:6]
                if (yyyymm_c >= sdt_yyyymm and yyyymm_c <= edt_yyyymm):
                    files_out.append(files[f])
    if (len(files_out) == 0):
        sys.exit('get_extract_files:  no files found!')
    return list(sorted(set(files_out)))

#---------------------------------------------------------------------------
# Get seasonal stats.
#---------------------------------------------------------------------------
def get_seasonal_files(files):

    springs = defaultdict(list)
    summers = defaultdict(list)
    falls   = defaultdict(list)
    winters = defaultdict(list)
    annuals = defaultdict(list)
    spring_ndays = {}
    summer_ndays = {}
    fall_ndays   = {}
    winter_ndays = {}
    annual_ndays = {}    
    years_all = []

    #--- Loop over all files, binning them into seasons / years.
    for f in range(len(files)):
        dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]
        if (dt_c):
            yyyy = dt_c[0:4]
            mm   = dt_c[4:6]
            years_all.append(yyyy)
            if (mm == '03' or mm == '04' or mm == '05'):
                if ((yyyy) in spring_ndays):
                    spring_ndays[yyyy] += 1
                else:
                    spring_ndays[yyyy] = 1
                springs[yyyy].append(files[f])
            elif (mm == '06' or mm == '07' or mm == '08'):
                if ((yyyy) in summer_ndays):
                    summer_ndays[yyyy] += 1
                else:
                    summer_ndays[yyyy] = 1
                summers[yyyy].append(files[f])
            elif (mm == '09' or mm == '10' or mm == '11'):
                if ((yyyy) in fall_ndays):
                    fall_ndays[yyyy] += 1
                else:
                    fall_ndays[yyyy] = 1
                falls[yyyy].append(files[f])
            elif (mm == '12'):
                yyyyw = str(int(yyyy) + 1)
                if ((yyyyw) in winter_ndays):
                    winter_ndays[yyyyw] += 1
                else:
                    winter_ndays[yyyyw] = 1
                winters[yyyyw].append(files[f])
            elif (mm == '01' or mm == '02'):
                if ((yyyy) in winter_ndays):
                    winter_ndays[yyyy] += 1
                else:
                    winter_ndays[yyyy] = 1
                winters[yyyy].append(files[f])

            if ((yyyy) in annual_ndays):
                annual_ndays[yyyy] += 1
            else:
                annual_ndays[yyyy] = 1
            annuals[yyyy].append(files[f])

    #--- Get list of unique, sorted years found above.
    years = list(sorted(set(years_all)))

    #--- Check number of data points per unique years and nan out
    #--- ones that do not have enough (< min_season_days).
    for y in range(len(years)):
        yyyy = years[y]
        if (yyyy in spring_ndays):
            if (spring_ndays[yyyy] < 3):
                springs[yyyy] = np.nan
        else:
            springs[yyyy] = np.nan
        if (yyyy in summer_ndays):
            if (summer_ndays[yyyy] < 3):
                summers[yyyy] = np.nan
        else:
            summers[yyyy] = np.nan
        if (yyyy in fall_ndays):            
            if (fall_ndays[yyyy] < 3):
                falls[yyyy] = np.nan
        else:
            falls[yyyy] = np.nan
        if (yyyy in winter_ndays):            
            if (winter_ndays[yyyy] < 3):
                winters[yyyy] = np.nan
        else:
            winters[yyyy] = np.nan
        if (yyyy in annual_ndays):            
            if (annual_ndays[yyyy] < 3):
                annuals[yyyy] = np.nan
        else:
            annuals[yyyy] = np.nan

    return springs, summers, falls, winters, annuals, years

#---------------------------------------------------------------------------
# Get seasonal files (new)
#---------------------------------------------------------------------------
def get_seasonal_files_new(files):

    seafiles = defaultdict(list)
    ndays = {}
    years_all = []

    #--- Loop over all files, binning them into seasons / years.
    for f in range(len(files)):
        dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]
        if (dt_c):
            yyyy = dt_c[0:4]
            mm   = dt_c[4:6]

            if mm == '12':
                yyyy = str(int(yyyy) + 1)
            else:
                yyyy = yyyy
            years_all.append(yyyy)

            if (mm == '03' or mm == '04' or mm == '05'):
                season = 'spring'
            if (mm == '06' or mm == '07' or mm == '08'):
                season = 'summer'
            if (mm == '09' or mm == '10' or mm == '11'):
                season = 'fall'
            if (mm == '12' or mm == '01' or mm == '02'):
                season = 'winter'

            if ((yyyy) in ndays):
                ndays[season+yyyy] += 1
            else:
                ndays[season+yyyy] = 1
            seafiles[season,yyyy].append(files[f])

    #--- Get list of unique, sorted years found above.
    years = list(sorted(set(years_all)))

    return seafiles, years

#---------------------------------------------------------------------------
# Loop over each models, grabbing data from extract wrfout files.
#---------------------------------------------------------------------------
def load_extract_data(geo_em, stns, latpts, lonpts, elevs, models, \
                      data_dirs, sdt, edt, varsin, pickle_dir, load_pickle):

    #--- Create temporary list of variables, removing 'T2MEAN' if in original
    #--- list (we'll calculate T2MEAN if necessary after loading data).
    vars_tmp = list(varsin)
    if 'T2MEAN' in vars_tmp:
        i = vars_tmp.index('T2MEAN')
        del vars_tmp[i]

    #--- Load geo em file to get lat/lon/elev grids.  Also interpolate
    #--- model elevation to stns lat/lon point locations.
    print 'Loading:\n', geo_em
    (lat, lon, hgt) = load_geo_em(geo_em)
    elevs_mod  = []
    elevs_diff = []    
    for ns in range(len(stns)):
        dat_interp = bilinear_interpolate(lat, lon, hgt, \
                                          latpts[ns], lonpts[ns], 0)
        elevs_mod.append(dat_interp)
        elevs_diff.append(elevs[ns] - dat_interp)

    for m in range(len(models)):
        data_all = {}
        dts_all = []
        model = models[m]
        pickle_file = get_mod_pkl_name(pickle_dir, model, sdt, edt)

        #--- Get list of extracted files for this model.
        print 'model = ', model
        files = get_extract_files(data_dirs, model, sdt, edt, 0)
        files = list(sorted(set(files)))

        print 'Loading: '
        #--- Loop over all extract files, reading them in.
        for f in range(len(files)):
    
            (vardat, ntimes, dt) = load_wrfout(files[f], vars_tmp)
            print files[f], ntimes
                
            #--- Loop over each vars_tmp, interpolating data to stns for
            #--- each ntimes.
            for v in range(len(vars_tmp)):
    
                var = vars_tmp[v]
                vardat_c = vardat[var]
                    
                #--- Loop over each ntime, grabbing interpolated data for
                #--- variable var for each stns.
                for n in range(ntimes):
                    dt_c = time_increment(dt, n, dtfmt)
                    for ns in range(len(stns)):
                        dat_interp = bilinear_interpolate(lat, lon, \
                                                          vardat_c[:,:,n],\
                                                          latpts[ns], \
                                                          lonpts[ns], 0)
    
                        #--- Do lapse rate correction for temperature.
                        if (var == 'T2MAX' or var == 'T2MIN'):
                            dat_c = dat_interp - std_lapse * elevs_diff[ns]
                        else:
                            dat_c = dat_interp
    
                        data_all[var, stns[ns], dt_c] = dat_c
                        dts_all.append(dt_c)
    
        #--- Get sorted unique list of all dates seen.
        dts_unique = list(sorted(set(dts_all)))

        #--- If T2MEAN was requested, calculate it from T2MAX and T2MIN for
        #--- every station and date.
        if 'T2MEAN' in varsin:
            for ns in range(len(stns)):
                for d in range(len(dts_unique)):
                    keymax = ('T2MAX', stns[ns], dts_unique[d])
                    keymin = ('T2MIN', stns[ns], dts_unique[d])
                    if keymax in data_all and keymin in data_all:
                        keymean = ('T2MEAN', stns[ns], dts_unique[d])

#                        if data_all[keymax] < (273-60) or \
#                           data_all[keymax] > (273+60) or \
#                           data_all[keymin] < (273-60) or \
#                           data_all[keymin] > (273+60):
#                            print 'stns[ns] = ', stns[ns]
#                            print 'dts_unique[d] = ', dts_unique[d]
#                            print 'max = ', data_all[keymax]
#                            print 'min = ', data_all[keymin]
##                            sys.exit()
                        
                        data_all[keymean] = np.mean([data_all[keymax], \
                                                     data_all[keymin]])

        #--- Create pickle file.
        pickle.dump((model, data_all, dts_unique, varsin, stns), \
                    open(pickle_file,'wb'), -1)

    return data_all, dts_unique, models, varsin, stns
