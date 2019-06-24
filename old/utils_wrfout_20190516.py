#!/usr/bin/python
import os, sys, glob, re, math
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from collections import defaultdict
from utils_cmap import *
import pickle
from utils_date import *
from utils_ghcnd_obs import *

dtfmt = "%Y%m%d%H"

#seasons    = ['winter']
#seasons    = ['annual', 'spring', 'summer', 'fall', 'winter']
seasons    = ['spring', 'summer', 'fall', 'winter', 'annual']
seasons_mo = ['ALL', 'MAM', 'JJA', 'SON', 'DJF']
seasons_lab = {
    'annual': 'Jan-Dec',    
    'spring': 'Mar-Apr-May',
    'summer': 'Jun-Jul-Aug',
    'fall':   'Sep-Oct-Nov',
    'winter': 'Dec-Jan-Feb'
    }
var_lab = {
    'T2MAX': 'Average Max Temperature ($^\circ$C)',
    'T2MIN': 'Average Min Temperature ($^\circ$C)',
    'PREC': 'Total Precipitation (in.)',
    'SPDUV10MEAN': 'Average Wind Speed (m/s)',
    'SPDUV10MAX': 'Maximum Wind Speed (m/s)',
    'SNOW': 'Snow Water Equivalent (mm)',
    'SWDOWN': 'Shortwave Down Radiation (W/m^2)'    
    }
ylims = {
    'PRECmin':  [0, 1.0],
    'PRECmax':  [0, 5.0],    
    'PRECavg':  [0, 45.0],
    'PRECtot':  [0, 40.0],        
    'T2MAXavg': [-5, 40.0],
    'T2MAXmax': [0, 50.0],    
    'T2MINmin': [-30.0, 25.0],
    'T2MINmax': [-30.0, 25.0],
    'T2MINavg': [-20, 25.0],
    'SPDUV10MEANavg': [0, 5.0],
    'SPDUV10MAXmax': [0, 20.0]    
    }
labels = {
    'PRECtot':  'Total Precipitation (in.)',
    'PRECmax':  'Maximum Precipitation (in.)',
    'PRECmin':  'Minimum Precipitation (in.)',
    'PRECavg':  'Average Precipitation (in.)',        
    'T2MAXavg': 'Average Maximum Temperature ($^\circ$C)',
    'T2MAXmax': 'Max Maximum Temperature ($^\circ$C)',
    'T2MAXmin': 'Min Maximum Temperature ($^\circ$C)',    
    'T2MINavg': 'Average Minimum Temperature ($^\circ$C)',
    'T2MINmax': 'Max Minimum Temperature ($^\circ$C)',
    'T2MINmin': 'Min Minimum Temperature ($^\circ$C)',
    'SPDUV10MEANavg': 'Average Wind Speed (m/s)',
    'SPDUV10MAXmax': 'Maximum Wind Speed (m/s)'    
    }

colorbar_labs = {
    'T2MAX': 'Temperature ($^\circ$C)',
    'T2MIN': 'Temperature ($^\circ$C)',
    'PREC': 'Precipitation (in.)',
    'SPDUV10MEAN': 'Wind Speed (m/s)',
    'SPDUV10MAX': 'Wind Speed (m/s)',
    'SNOW': 'Snow Water Equivalent (mm)',
    'SWDOWN': 'Shortwave Down Radiation (W/m^2)'     
    }

zooms = ['Z1', 'Z2']
target_lat = 47.44472
target_lon = -122.31361
fudge_lat = 1.7
fudge_lon = 2.1
minlatZ2 = target_lat - fudge_lat
maxlatZ2 = target_lat + fudge_lat
minlonZ2 = target_lon - fudge_lon
maxlonZ2 = target_lon + fudge_lon
minlat = {
    'Z1': 39.0,
    'Z2': minlatZ2 }
maxlat = {
    'Z1': 51.0,
    'Z2': maxlatZ2 }
minlon = {
    'Z1': -131.0,
    'Z2': minlonZ2 }
maxlon = {
    'Z1': -108.0,
    'Z2': maxlonZ2 }

mm2in           = 0.0393700787
min_season_days = 85
std_lapse       = 0.0065  #-- std. atmos. lapse rate
smooth_fact     = 5       #-- odd number, used for smoothing time series plots.

fs      = 9
titlefs = 9
width   = 10
height  = 8
maplw   = 1.0
mvc     = -9999.0

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
# Load wrfout file.
#---------------------------------------------------------------------------
def sum_wrfout_slice(files, vars, maski, maskj):    
    sums = {}
    cnts = {}

    for f in range(len(files)):
        print 'reading ', files[f]
        nc = NetCDFFile(files[f], 'r')
        ntimes = len(nc.dimensions['Time'])

        for v in range(len(vars)):
            var = vars[v]
            dat_tmp = nc.variables[var][:,:,:]
            dat_tmp = np.transpose(dat_tmp)

            sum_c = 0
            cnt_c = 0
            for n in range(len(maski)):
                sum_c = sum_c + sum(dat_tmp[maski[n], maskj[n],:])
                cnt_c = cnt_c + len(dat_tmp[maski[n], maskj[n],:])

            if (var in sums):
                sums[var] = sums[var] + sum_c
                cnts[var] = cnts[var] + cnt_c
            else:
                sums[var] = sum_c
                cnts[var] = cnt_c

    return sums, cnts

#---------------------------------------------------------------------------
# For ensemble standard deviation plots....
#---------------------------------------------------------------------------
def get_ensmean_std(file_plot, var, models, years, season, files):

    seafiles = defaultdict(list)
    ndays = {}

#    springs = defaultdict(list)
#    summers = defaultdict(list)
#    falls   = defaultdict(list)
#    winters = defaultdict(list)
#    spring_ndays = {}
#    summer_ndays = {}
#    fall_ndays   = {}
#    winter_ndays = {}

    years_all = []
    dts_all = []
    #--- Loop over all files, binning them into seasons / years.
    for f in range(len(files)):
        
        dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]
        print files[f], dt_c
        
#        if (dt_c):
#            dts_all.append(dt_c)
#            yyyy = dt_c[0:4]
#            mm   = dt_c[4:6]
#
#            if mm == '12':
#                yyyy = str(int(yyyy) + 1)
#            else:
#                yyyy = yyyy
#            years_all.append(yyyy)
#
#            if (mm == '03' or mm == '04' or mm == '05'):
#                season = 'spring'
#            if (mm == '06' or mm == '07' or mm == '08'):
#                season = 'summer'
#            if (mm == '09' or mm == '10' or mm == '11'):
#                season = 'fall'
#            if (mm == '12' or mm == '01' or mm == '02'):
#                season = 'winter'
#
#            if ((yyyy) in ndays):
#                ndays[season+yyyy] += 1
#            else:
#                ndays[season+yyyy] = 1
#            seafiles[season,yyyy].append(files[f])

    sys.exit()

    
#    syscom = ncocom
#    nyears = 0
#    nmodels = len(models)
#    files = []
#    dts_all = []
#    for y in range(len(years)):
#        yyyy = years[y]
#
#        for m in range(nmodels):
#            model = models[m]
#            sum_files_c = sum_files[model]
#            spr_files_c = spr_files[model]
#            win_files_c = win_files[model]
#            fal_files_c = fal_files[model]
#            ann_files_c = ann_files[model]
#
#            if season == 'summer':
#                if (sum_files_c[yyyy] == sum_files_c[yyyy]):
#                    nyears += 1
#                    for f in range(len(sum_files_c[yyyy])):
#                        files.append(sum_files_c[yyyy][f])
#                    else:
#                        continue
#            elif season == 'spring':
#                if (spr_files_c[yyyy] == spr_files_c[yyyy]):
#                    nyears += 1
#                    for f in range(len(spr_files_c[yyyy])):
#                        files.append(spr_files_c[yyyy][f])
#                    else:
#                        continue
#            elif season == 'winter':
#                if (win_files_c[yyyy] == win_files_c[yyyy]):
#                    nyears += 1
#                    for f in range(len(win_files_c[yyyy])):
#                        files.append(win_files_c[yyyy][f])
#                    else:
#                        continue
#            elif season == 'fall':                        
#                if (fal_files_c[yyyy] == fal_files_c[yyyy]):
#                    nyears += 1
#                    for f in range(len(fal_files_c[yyyy])):
#                        files.append(fal_files_c[yyyy][f])
#                    else:
#                        continue
#            elif season == 'annual':
#                if (ann_files_c[yyyy] == ann_files_c[yyyy]):
#                    nyears += 1
#                    for f in range(len(ann_files_c[yyyy])):
#                        files.append(ann_files_c[yyyy][f])
#                    else:
#                        continue
#
#    #--- Read in files from above.
#    dat_all = np.ones((nx, ny, ndts, nmodels)) * np.nan
#    
##    nfiles = len(files)
#    nfiles = 10
#    for f in range(nfiles):
#        print files[f]
#        (vardat, ntimes, dt) = load_wrfout(files[f], [var])
#        plotdat = vardat[var][:,:,0]
#        nx,ny = plotdat.shape
#        if f == 0:
#            dat_all = np.ones((nx, ny, nfiles)) * np.nan
#        dat_all[:,:,f] = plotdat
#
#    test = np.std(dat_all, axis=2)
#    print test.shape
#    print test
#    sys.exit()

    return files, nyears

#---------------------------------------------------------------------------
# For ensemble mean plots, get system nco command.
#---------------------------------------------------------------------------
def get_ensmean_syscom(file_plot, var, models, years, season, \
                       sum_files, spr_files, \
                       win_files, fal_files, ann_files, ncocom):
    syscom = ncocom
    nyears = 0
    for m in range(len(models)):
        model = models[m]
        sum_files_c = sum_files[model]
        spr_files_c = spr_files[model]
        win_files_c = win_files[model]
        fal_files_c = fal_files[model]
        ann_files_c = ann_files[model]
        for y in range(len(years)):
            yyyy = years[y]
            if season == 'summer':
                if (sum_files_c[yyyy] == sum_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(sum_files_c[yyyy])):
                        syscom = syscom + ' ' + sum_files_c[yyyy][f]
                    else:
                        continue
            elif season == 'spring':
                if (spr_files_c[yyyy] == spr_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(spr_files_c[yyyy])):
                        syscom = syscom + ' ' + spr_files_c[yyyy][f]
                    else:
                        continue
            elif season == 'winter':
                if (win_files_c[yyyy] == win_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(win_files_c[yyyy])):
                        syscom = syscom + ' ' + win_files_c[yyyy][f]
                    else:
                        continue
            elif season == 'fall':                        
                if (fal_files_c[yyyy] == fal_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(fal_files_c[yyyy])):
                        syscom = syscom + ' ' + fal_files_c[yyyy][f]
                    else:
                        continue
            elif season == 'annual':
                if (ann_files_c[yyyy] == ann_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(ann_files_c[yyyy])):
                        syscom = syscom + ' ' + ann_files_c[yyyy][f]
                    else:
                        continue

    return syscom, nyears

#---------------------------------------------------------------------------
# For ensemble mean plots, run system nco command.
#---------------------------------------------------------------------------
def run_ensmean_syscom(file_plot, var, ncocom, syscom, yyyy_s, yyyy_e, \
                       season, rundir, ):

    #--- If system command is too big (too many files to send
    #--- to nco), break it up into 1000 piece chunks.
    files_c = syscom.replace(ncocom, '')
    scsplit = files_c.split()
    nfiles = len(scsplit)

    niter = 1000
    if nfiles > niter:
        nt =  int(math.ceil(float(nfiles) / float(niter)))
        files_all = []
        for n in range(nt):
            file_c = rundir + '/' + var + '_' + season + '_' + \
                     yyyy_s + '_' + yyyy_e + '_subset' + str(n) + '.nc'
            if os.path.isfile(file_c):
                print 'this file already exists: ', file_c
                files_all.append(file_c)
                continue
            nb = n * niter
            ne = (n+1) * niter
            syscom = ncocom + ' ' + ' '.join(scsplit[nb:ne]) + ' ' + file_c
            print syscom
            iret = os.system(syscom)
            if iret != 0:
                print 'problem running nco! iret = ', iret
                return mvc
            files_all.append(file_c)
        syscom = ncocom + ' ' + ' '.join(files_all) + ' ' + file_plot
        print syscom                
        iret = os.system(syscom)
        if iret != 0:
            print 'problem running nco! iret = ', iret
            return mvc
    else:
        syscom = syscom + ' ' + file_plot
        print 'syscom = ', syscom
        iret = os.system(syscom)
        if iret != 0:
            print 'problem running nco! iret = ', iret
            return mvc

    return 1

#---------------------------------------------------------------------------
# Get plot file name.
#---------------------------------------------------------------------------
def get_plotfname(plot_dir, stn, sdt, edt, season, var, stat):
    plotfname = plot_dir + '/' + stn + '_' + sdt + '_' + edt + '_' + \
                season + '_' + var + '_' + stat + '.png'
    return plotfname

#---------------------------------------------------------------------------
# Get plot title.
#---------------------------------------------------------------------------
def get_title(season, stn, var, stat, years):
    titleout = station_name_dict[stn] + ' (' + stn.upper() + '), ' + \
               seasons_lab[season] + ' ' + labels[var+stat] + ', ' + \
               str(years[0]) + ' - ' + str(years[len(years)-1])
    return titleout

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

#    springs = defaultdict(list)
#    summers = defaultdict(list)
#    falls   = defaultdict(list)
#    winters = defaultdict(list)
#    spring_ndays = {}
#    summer_ndays = {}
#    fall_ndays   = {}
#    winter_ndays = {}

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

'''
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
'''


#---------------------------------------------------------------------------
# Get seasonal stats (totals or averages currently).
#---------------------------------------------------------------------------
def get_seasonal_stats(data_all, dts_all, models, stns, var, stat):

    out_sum = {}
    out_avg = {}
    out_max = {}
    out_min = {}
        
    years_all = []
    stn_pr = 'KSEA'
    
    for m in range(len(models)):
        mod = models[m]
        data_all_c = data_all[mod]
        for s in range(len(stns)):
            stn = stns[s]
            ndays = {}
            years_mod_stn = []
            for d in range(len(dts_all)):
                yyyy = dts_all[d][0:4]
                mm = dts_all[d][4:6]
                years_mod_stn.append(yyyy)

                if ((var,stn,dts_all[d]) in data_all_c):
                    dat_c = data_all_c[var,stn,dts_all[d]]
                else:
                    ndays[season+yyyy] = np.nan
                    out_sum[season,mod,stn,yyyy] = np.nan
                    out_max[season,mod,stn,yyyy] = np.nan
                    out_min[season,mod,stn,yyyy] = np.nan
                    continue
                
                #--- Spring.
                if (mm == '03' or mm == '04' or mm == '05'):
                    season = 'spring'
                    #--- Count number of days for this year's spring.
                    if (season+yyyy in ndays):
                        ndays[season+yyyy] += 1
                    else:
                        ndays[season+yyyy] = 1

                    if ((season,mod,stn,yyyy) in out_sum):
                        out_sum[season,mod,stn,yyyy] = out_sum[\
                            season,mod,stn,yyyy] + dat_c
                    else:
                        out_sum[season,mod,stn,yyyy] = dat_c

                    #--- Get maximum.
                    if ((season,mod,stn,yyyy) in out_max):
                        if (dat_c > out_max[season,mod,stn,yyyy]):
                            out_max[season,mod,stn,yyyy] = dat_c

#                            if (stn == stn_pr):
#                                print '-----'
#                                print stn, mod, yyyy, dts_all[d], dat_c
#                                print 'out_max[season,mod,stn,yyyy] = ', \
#                                      out_max[season,mod,stn,yyyy]
                    else:
                        out_max[season,mod,stn,yyyy] = dat_c

                    #--- Get minimum.
                    if ((season,mod,stn,yyyy) in out_min):
                        if (dat_c < out_min[season,mod,stn,yyyy]):
                            out_min[season,mod,stn,yyyy] = dat_c
                    else:
                        out_min[season,mod,stn,yyyy] = dat_c
                        
                #--- Summer.
                if (mm == '06' or mm == '07' or mm == '08'):
                    season = 'summer'
                    #--- Count number of days for this year's summer.
                    if (season+yyyy in ndays):
                        ndays[season+yyyy] += 1
                    else:
                        ndays[season+yyyy] = 1
                    if ((season,mod,stn,yyyy) in out_sum):
                        out_sum[season,mod,stn,yyyy] = out_sum[\
                            season,mod,stn,yyyy] + dat_c
                    else:
                        out_sum[season,mod,stn,yyyy] = dat_c

                    #--- Get maximum.
                    if ((season,mod,stn,yyyy) in out_max):
                        if (dat_c > out_max[season,mod,stn,yyyy]):
                            out_max[season,mod,stn,yyyy] = dat_c
                    else:
                        out_max[season,mod,stn,yyyy] = dat_c

                    #--- Get minimum.
                    if ((season,mod,stn,yyyy) in out_min):
                        if (dat_c < out_min[season,mod,stn,yyyy]):
                            out_min[season,mod,stn,yyyy] = dat_c
                    else:
                        out_min[season,mod,stn,yyyy] = dat_c

                #--- Fall.
                if (mm == '09' or mm == '10' or mm == '11'):
                    season = 'fall'
                    #--- Count number of days for this year's fall.
                    if (season+yyyy in ndays):
                        ndays[season+yyyy] += 1
                    else:
                        ndays[season+yyyy] = 1
                    if ((season,mod,stn,yyyy) in out_sum):
                        out_sum[season,mod,stn,yyyy] = out_sum[\
                            season,mod,stn,yyyy] + dat_c
                    else:
                        out_sum[season,mod,stn,yyyy] = dat_c

                    #--- Get maximum.
                    if ((season,mod,stn,yyyy) in out_max):
                        if (dat_c > out_max[season,mod,stn,yyyy]):
                            out_max[season,mod,stn,yyyy] = dat_c
                    else:
                        out_max[season,mod,stn,yyyy] = dat_c

                    #--- Get minimum.
                    if ((season,mod,stn,yyyy) in out_min):
                        if (dat_c < out_min[season,mod,stn,yyyy]):
                            out_min[season,mod,stn,yyyy] = dat_c
                    else:
                        out_min[season,mod,stn,yyyy] = dat_c

                #--- Winter, December only.  Put December winter data into
                #--- next year's winter-- naming winters by their Jan/Feb
                #--- year rather than their Dec year.
                if (mm == '12'):
                    season = 'winter'
                    yyyyw = str(int(yyyy) + 1)
                    #--- Count number of days for this year's winter.
                    if (season+yyyyw in ndays):
                        ndays[season+yyyyw] += 1
                    else:
                        ndays[season+yyyyw] = 1
                    if ((season,mod,stn,yyyyw) in out_sum):
                        out_sum[season,mod,stn,yyyyw] = out_sum[\
                            season,mod,stn,yyyyw] + dat_c
                    else:
                        out_sum[season,mod,stn,yyyyw] = dat_c

                    #--- Get maximum.
                    if ((season,mod,stn,yyyyw) in out_max):
                        if (dat_c > out_max[season,mod,stn,yyyyw]):
                            out_max[season,mod,stn,yyyyw] = dat_c
                    else:
                        out_max[season,mod,stn,yyyyw] = dat_c

                    #--- Get minimum.
                    if ((season,mod,stn,yyyyw) in out_min):
                        if (dat_c < out_min[season,mod,stn,yyyyw]):
                            out_min[season,mod,stn,yyyyw] = dat_c
                    else:
                        out_min[season,mod,stn,yyyyw] = dat_c
                        
                #--- Winter, Jan & Feb.
                if (mm == '01' or mm == '02'):
                    season = 'winter'
                    #--- Count number of days for this year's winter.
                    if (season+yyyy in ndays):
                        ndays[season+yyyy] += 1
                    else:
                        ndays[season+yyyy] = 1
                    if ((season,mod,stn,yyyy) in out_sum):
                        out_sum[season,mod,stn,yyyy] = out_sum[\
                            season,mod,stn,yyyy] + dat_c
                    else:
                        out_sum[season,mod,stn,yyyy] = dat_c

                    #--- Get maximum.
                    if ((season,mod,stn,yyyy) in out_max):
                        if (dat_c > out_max[season,mod,stn,yyyy]):
                            out_max[season,mod,stn,yyyy] = dat_c
                    else:
                        out_max[season,mod,stn,yyyy] = dat_c

                    #--- Get minimum.
                    if ((season,mod,stn,yyyy) in out_min):
                        if (dat_c < out_min[season,mod,stn,yyyy]):
                            out_min[season,mod,stn,yyyy] = dat_c
                    else:
                        out_min[season,mod,stn,yyyy] = dat_c
                        
            #--- Get list of unique, sorted years found above.
            years = list(sorted(set(years_mod_stn)))

#            if (mod == 'bcc-csm1.1' and stn == 'KSEA'):
#                for yyyy in years:
#                    if ((season+yyyy) in ndays):
#                        print season,yyyy, ndays[season+yyyy]
#                    else:
#                        print season,yyyy, 'missing'

            #--- Check number of data points per unique years and nan out
            #--- ones that do not have enough (< min_season_days).
            for y in range(len(years)):
                yyyy = years[y]
                for s in range(len(seasons)):
                    season = seasons[s]
                    if ((season+yyyy) in ndays):
                        if (ndays[season+yyyy] < min_season_days):
                            out_sum[season,mod,stn,yyyy] = np.nan
                            out_max[season,mod,stn,yyyy] = np.nan
                            out_min[season,mod,stn,yyyy] = np.nan
                    else:
                        ndays[season+yyyy] = np.nan
                        out_sum[season,mod,stn,yyyy] = np.nan
                        out_max[season,mod,stn,yyyy] = np.nan
                        out_min[season,mod,stn,yyyy] = np.nan

            #--- if stat requested was 'avg', do averaging.
            if (stat == 'avg' or stat == 'max' or stat == 'min'):
                for y in range(len(years)):
                    yyyy = years[y]
                    for s in range(len(seasons)):
                        season = seasons[s]
                        out_avg[season,mod,stn,yyyy] = out_sum[\
                            season,mod,stn,yyyy] / ndays[season+yyyy]
                
            #--- Keep around all years found for this model and station.
            years_all.append(years_mod_stn)

    #--- Final, unique, sorted list of years found.
    years = list(sorted(set(years_all[0])))            

    if (stat == 'avg'):
        return out_avg, years
    elif (stat == 'max'):
        return out_max, years
    elif (stat == 'min'):
        return out_min, years
    else:
        return out_sum, years

#---------------------------------------------------------------------------
# Make time series of sent-in data for a set of stations for one model.
#---------------------------------------------------------------------------
def ts_stns(statplot, years, stns, mod, titlein, plotfname, var, stat, \
            cols, seasonsplot):

    fig, ax = plt.subplots( figsize=(10,6) )
    xlabs = []
    for s in range(len(stns)):
        ptot = []
        stn = stns[s]
        for y in range(len(years)):
            yyyy = years[y]
            for ss in range(len(seasons)):
                season = seasons[ss]
                if (season in seasonsplot):
                    if (var == 'PREC'):
                        ptot.append(statplot[season,mod,stns[s],yyyy] * mm2in)
                    elif (var == 'T2MAX' or var == 'T2MIN'):
                        ptot.append(statplot[season,mod,stns[s],yyyy] - 273.15)
                    elif (var.find('SPD') >= 0):
                        ptot.append(statplot[season,mod,stns[s],yyyy])
                    if (s == 0):
                        xlabs.append(season.title() + ' ' + years[y])
        plt.plot(ptot, alpha = 0.3, label = stn, color = cols[s])
        lab_c = stn + ', ' + str(smooth_fact) + '-pt smoothing'
        plt.plot(smooth(ptot, smooth_fact), label = lab_c, color = cols[s])

    #--- y-axis labeling.
    plt.ylim(ylims[var+stat])
    plt.ylabel('Seasonal ' + labels[var+stat], fontsize=fs+1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    

    #--- x-axis labels.
    xticks_c = range(0,len(xlabs), 4)
    xlabs_c = [ xlabs[i] for i in xticks_c ]    
    plt.xticks(xticks_c)
    plt.tick_params(axis='x', which='major', labelsize=fs-1)
    ax.set_xticklabels(xlabs_c, rotation=90)        
    plt.xlabel('Season', fontsize=fs+1)

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    plt.tight_layout()
    plt.grid()
    plt.legend(fontsize = fs-1, loc = 'best')

    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)
    plt.close()

    return 1

#---------------------------------------------------------------------------
# Make time series of sent-in data for a set of models for one station.
#---------------------------------------------------------------------------
def ts_mods(statplot, years, obsplot, years_obs, stn, models, titlein, \
            plotfname, var, stat, cols, seasonsplot):

    fig, ax = plt.subplots( figsize=(10,6) )

    lw_sm = 1.8

    #--- Plot models.
    xlabs = []
    mod_all = np.ones((len(models),len(years))) * np.nan
    for m in range(len(models)):
        ptot = []
        mod = models[m]
        for y in range(len(years)):
            yyyy = years[y]
            for ss in range(len(seasons)):
                season = seasons[ss]
                key = (season,mod,stn,yyyy)
                if (season in seasonsplot):
                    if (var == 'PREC'):
                        ptot_c = statplot[key] * mm2in
                    elif (var == 'T2MAX' or var == 'T2MIN'):
                        ptot_c = statplot[key] - 273.15
                    elif (var.find('SPD') >= 0):
                        ptot_c = statplot[key]
                    if (m == 0):
                        xlabs.append(years[y])

                    ptot.append(ptot_c)
                    
            mod_all[m,y] = ptot_c
            
#        plt.plot(ptot, alpha = 0.20, label = mod.upper(), color = cols[m], \
#                 linewidth=lw_sm)
        lab_c = mod.upper() + ', ' + str(smooth_fact) + '-pt smoothing'
        plt.plot(smooth(ptot, smooth_fact), label = lab_c, color = cols[m], \
                 linewidth=lw_sm, alpha = 0.23)

    #--- Plot ensemble mean.
    plt.plot(np.mean(mod_all, axis=0), label = 'Ensemble Mean', \
             color = 'darkgreen', linewidth=2.2)
    
    #--- Plot observations.
    otot = []
    for y in range(len(years_obs)):
#        yyyy = years[y]
        yyyy = years_obs[y]
        for ss in range(len(seasons)):
            season = seasons[ss]
            key = (season,stn,yyyy)
            if (season in seasonsplot):
                otot.append(obsplot[key])
    plt.plot(otot, alpha=0.5, label='Observed', color='black', \
             linestyle='None', marker='o', markersize=3)
    lab_c = 'Observed ' + str(smooth_fact) + '-pt smoothing'
    plt.plot(smooth(otot,smooth_fact), label=lab_c, color='black', \
             linewidth=lw_sm)

    #--- y-axis labeling.
    plt.ylim(ylims[var+stat])
    plt.ylabel('Seasonal ' + labels[var+stat], fontsize=fs+1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    

    #--- x-axis labels.
    xticks_c = range(0,len(xlabs), 4)
    xlabs_c = [ xlabs[i] for i in xticks_c ]    
    plt.xticks(xticks_c)
    plt.tick_params(axis='x', which='major', labelsize=fs-1)
    ax.set_xticklabels(xlabs_c, rotation=90)        
    plt.xlabel('Year', fontsize=fs+1)

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    plt.tight_layout()
    plt.grid()
    plt.legend(fontsize = fs-1, loc = 'best')

    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)
    plt.close()

    return 1

#---------------------------------------------------------------------------
# Loop over each models, grabbing data from extract wrfout files.
#---------------------------------------------------------------------------
def load_extract_data(geo_em, stns, latpts, lonpts, elevs, models, \
                      data_dirs, sdt, edt, vars, pickle_dir, load_pickle):

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
        pickle_file = pickle_dir + '/extract_' + model + '_' + \
                      sdt + '_' + edt + '.pickle'

        #--- Get list of extracted files for this model.
        files = get_extract_files(data_dirs, model, sdt, edt, 0)

        files = list(sorted(set(files)))

        print 'Loading: '
        #--- Loop over all extract files, reading them in.
        for f in range(len(files)):
    
            (vardat, ntimes, dt) = load_wrfout(files[f], vars)
            print files[f], ntimes
                
            #--- Loop over each vars, interpolating data to stns for
            #--- each ntimes.
            for v in range(len(vars)):
    
                var = vars[v]
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
        
        #--- Create pickle file.
        pickle.dump((model, data_all, dts_unique, vars, stns), \
                    open(pickle_file,'wb'), -1)

    return data_all, dts_unique, models, vars, stns

#---------------------------------------------------------------------------
# Interpolate from a sent-in grid to a sent in lat/lon point.
#---------------------------------------------------------------------------
def bilinear_interpolate(latgrid, longrid, datagrid, latpt, lonpt, deb):

    totdiff = abs(latgrid - latpt) + abs(longrid - lonpt)
    (i,j) = np.unravel_index(totdiff.argmin(), totdiff.shape)

    #--- Get lat/lon box in which latpt,lonpt resides.
    if (latpt >= latgrid[i,j] and lonpt >= longrid[i,j]):
        iif = i
        jjf = j
    elif (latpt < latgrid[i,j] and lonpt < longrid[i,j]):
        iif = i-1
        jjf = j-1
    elif (latpt >= latgrid[i,j] and lonpt < longrid[i,j]):
        iif = i-1
        jjf = j
    elif (latpt < latgrid[i,j] and lonpt >= longrid[i,j]):
        iif = i
        jjf = j-1

    (nx, ny) = np.shape(latgrid)

    if (deb == 1):
        print 'nx, ny  = ', nx, ny
        print 'iif,jjf = ', iif,jjf
        print 'latgrid[iif+1,jjf] = ', latgrid[iif+1,jjf]
        print 'latgrid[iif,jjf]   = ', latgrid[iif,jjf]

    if (iif >= (nx-1) or jjf >= (ny-1) or iif < 0 or jjf < 0):
        return mvc

    #--- Do bilinear interpolation to latpt,lonpt.
    dlat  = latgrid[iif,jjf+1] - latgrid[iif,jjf]
    dlon  = longrid[iif+1,jjf] - longrid[iif,jjf]
    dslat = latgrid[iif,jjf+1] - latpt
    dslon = longrid[iif+1,jjf] - lonpt

    wrgt = 1 - (dslon/dlon)
    wup  = 1 - (dslat/dlat)

    vll = datagrid[iif,jjf]
    vlr = datagrid[iif,jjf+1]
    vul = datagrid[iif+1,jjf]
    vur = datagrid[iif+1,jjf+1]

    if (deb > 1):
        print 'll lat, lon, val = ', latgrid[iif,jjf], longrid[iif,jjf], vll
        print 'lr lat, lon, val = ', latgrid[iif+1,jjf], longrid[iif+1,jjf], vlr
        print 'ur lat, lon, val = ', latgrid[iif+1,jjf+1], longrid[iif+1,jjf+1],vur
        print 'ul lat, lon, val = ', latgrid[iif,jjf+1], longrid[iif,jjf+1], vul
        print 'latpt, lonpt = ', latpt, lonpt
        
        print 'vll = ', vll
        print 'vlr = ', vlr
        print 'vul = ', vul
        print 'vur = ', vur

    datout = (1-wrgt) * ((1-wup) * vll + wup * vul) + \
        (wrgt) * ((1-wup) * vlr + wup * vur)

#    if (deb == 1):
#        print 'datout = ', datout
#        sys.exit()
#    if (datout == 0.0):
#        print 'nx, ny  = ', nx, ny
#        print 'iif,jjf = ', iif,jjf
#        print 'latgrid[iif+1,jjf] = ', latgrid[iif+1,jjf]
#        print 'latgrid[iif,jjf]   = ', latgrid[iif,jjf]
#        
#        print 'll lat, lon, val = ', latgrid[iif,jjf], longrid[iif,jjf], vll
#        print 'lr lat, lon, val = ', latgrid[iif+1,jjf], longrid[iif+1,jjf], vl#r
#        print 'ur lat, lon, val = ', latgrid[iif+1,jjf+1], longrid[iif+1,jjf+1],vur
#        print 'ul lat, lon, val = ', latgrid[iif,jjf+1], longrid[iif,jjf+1], vu#l
#        print 'latpt, lonpt = ', latpt, lonpt
#        
#        print 'vll = ', vll
#        print 'vlr = ', vlr
#        print 'vul = ', vul
#        print 'vur = ', vur
#
#        sys.exit()

    return datout

#---------------------------------------------------------------------------
# Run mapper.
#---------------------------------------------------------------------------
def run_mapper(var, model, rundir, season, files, yyyy, \
               lat, lon, levs_c, my_cmap ):

    file_plot = rundir + '/' + season + '_' + yyyy + '_' + var + '.nc'
    if (var == 'T2MAX' or var == 'T2MIN'):
        if (os.path.exists(file_plot)):
            print 'nc file ', file_plot, ' already exists-- using it'
        else:
            syscom = 'ncra'
            for f in range(len(files[yyyy])):
                print yyyy, files[yyyy][f]
                syscom = syscom + ' ' + files[yyyy][f]
            syscom = syscom + ' ' + file_plot
            print 'syscom = ', syscom        
            os.system(syscom)
    elif (var == 'PREC'):
        if (os.path.exists(file_plot)):
            print 'nc file ', file_plot, ' already exists-- using it'
        else:
            syscom = 'ncra -y ttl'
            for f in range(len(files[yyyy])):
                print yyyy, files[yyyy][f]
                syscom = syscom + ' ' + files[yyyy][f]
            syscom = syscom + ' ' + file_plot
            print 'syscom = ', syscom        
            os.system(syscom)

    print 'Loading ', var, ' from ', file_plot
    (vardat, ntimes, dt) = load_wrfout(file_plot, [var])

    titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
              seasons_lab[season] + ' ' + yyyy

    plotdat = vardat[var][:,:,0]
    if (var == 'PREC'):
        plotdat = plotdat * mm2in

    for z in range(len(zooms)):
        plotfname = rundir + '/' + model + '_' + yyyy + '_' + season + '_' + \
                    var + '_' + zooms[z] + '.png'
        print 'colorbar_labs[var] = ', colorbar_labs[var]
        mapper(var, lat, lon, plotdat, levs_c, my_cmap, maplw, \
               colorbar_labs[var], titlein, plotfname, zooms[z])
    return 1

#---------------------------------------------------------------------------
# Mapper.
#---------------------------------------------------------------------------
def mapper(var, lat, lon, grid, levs, cmap_in, maplw, colorbar_lab, titlein, \
           plotfname, zoom):

    ur_lat = maxlat[zoom]
    ur_lon = maxlon[zoom]
    ll_lat = minlat[zoom]
    ll_lon = minlon[zoom]
    lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
    lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

    if (zoom == 'Z1'):
        res = 'i'
    elif (zoom == 'Z2'):
        res = 'h'
        
    fig = plt.figure(figsize=(width,height))
    # left, bottom, width, height:
    ax = fig.add_axes([0.00,0.05,0.99,0.91])
    map = Basemap(resolution = res,projection='lcc',\
                  llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
                  urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
                  lat_0=lat_ctr,lon_0=lon_ctr,lat_1=(ur_lat - ll_lat))

    #--- Get lat and lon data in map's x/y coordinates.
    x,y = map(lon, lat)

    #--- Draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth = maplw)
    map.drawstates(linewidth = maplw)
    map.drawcountries(linewidth = maplw)

    #--- Draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(linewidth = maplw)

    #--- Draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0, 360, 30), linewidth = maplw)
    map.drawparallels(np.arange(-90, 90, 30), linewidth = maplw)

    print 'cmap = ', cmap_in
    print 'levs = ', levs
#    print 'cmaps[var] = ', cmaps[var]
#    print 'cmap_temp = ', cmap_temp

    print 'x.shape = ', x.shape
    print 'grid.shape = ', grid.shape    

    cs = plt.contourf(x, y, grid, levs, cmap=cmap_in)
#    cs = plt.contourf(x, y, grid)

    cbar = map.colorbar(cs, location='bottom', pad="3%", size=0.1, ticks=levs)
    cbar.set_label(colorbar_lab, fontsize=fs, size=fs-1)
    cbar.ax.tick_params(labelsize=fs-2)

    csl = plt.contour(x, y, grid, levs, colors = 'black', linewidths=0.8)
#    csl = plt.contour(x, y, grid, colors = 'black', linewidths=0.8)    
    plt.clabel(csl, inline=1, fontsize=10, fmt='%2.1f', fontweights='bold')

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    #--- Save plot.
    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)

    plt.close()

#---------------------------------------------------------------------------
# Gaussian smoother, taken from
# https://www.swharden.com/wp/2008-11-17-linear-data-smoothing-in-python/
#---------------------------------------------------------------------------
def smoothListGaussian(list,strippedXs=False,degree=5):  

     window=degree*2-1  
     weight=np.array([1.0]*window)  
     weightGauss=[]  

     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(np.exp((4*(frac))**2))  
         weightGauss.append(gauss)  

     weight=np.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  

     for i in range(len(smoothed)):  
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  

     return smoothed  

#def running_mean(x, N):
#    cumsum = np.cumsum(np.insert(x, 0, 0)) 
#    return (cumsum[N:] - cumsum[:-N]) / float(N)

#---------------------------------------------------------------------------
# Smoother like Matlab's running average smoother 'smooth', taken from
# https://stackoverflow.com/questions/40443020/matlabs-smooth-implementation-n-point-moving-average-in-numpy-python
#---------------------------------------------------------------------------
def smooth(a,WSZ):
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))
