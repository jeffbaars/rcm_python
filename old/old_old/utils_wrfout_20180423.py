#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from collections import defaultdict
from utils_cmap import *

seasons = ['spring', 'summer', 'fall', 'winter']
seasons_lab = {
    'spring': 'Spring (MAM)',
    'summer': 'Summer (JJA)',
    'fall':   'Fall (SON)',
    'winter': 'Winter (DJF)'
    }
var_lab = {
    'T2MAX': 'Average Max Temperature (K)',
    'T2MIN': 'Average Min Temperature (K)',
    'PREC': 'Total Precipitation (in.)',
    'SPDUV10MEAN': 'Average Wind Speed (m/s)',
    'SPDUV10MAX': 'Maximum Wind Speed (m/s)'    
    }
ylims = {
    'PRECmax':  [0, 5.0],
    'PRECtot':  [0, 40.0],    
    'T2MAXavg': [-5, 40.0],
    'T2MAXmax': [0, 50.0],    
    'T2MINmin': [-30.0, 25.0],
    'T2MINavg': [-20, 25.0],
    'SPDUV10MEANavg': [0, 5.0],
    'SPDUV10MAXmax': [0, 20.0]    
    }
ylabs = {
    'PRECtot': 'Seasonal Total Precipitation (in.)',
    'PRECmax': 'Seasonal Maximum Precipitation (in.)',
    'T2MAXavg': 'Seasonal Average Maximum Temperature ($^\circ$C)',
    'T2MAXmax': 'Seasonal Max Maximum Temperature ($^\circ$C)',
    'T2MINavg': 'Seasonal Average Minimum Temperature ($^\circ$C)',
    'T2MINmax': 'Seasonal Max Minimum Temperature ($^\circ$C)',
    'T2MINmin': 'Seasonal Min Minimum Temperature ($^\circ$C)',
    'SPDUV10MEANavg': 'Seasonal Average Wind Speed (m/s)',
    'SPDUV10MAXmax': 'Seasonal Maximum Wind Speed (m/s)'    
    }
colorbar_labs = {
    'T2MAX': 'Temperature (K)',
    'T2MIN': 'Temperature (K)',
    'PREC': 'Precipitation (in.)',
    'SPDUV10MEAN': 'Wind Speed (m/s)',
    'SPDUV10MAX': 'Wind Speed (m/s)'    
    }

zooms = ['Z1', 'Z2']
target_lat = 47.44472
target_lon = -122.31361
fudge_lat = 1.5
fudge_lon = 2.0
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
# Get list of extracted wrfout files between sdt_yyyymm and edt_yyyymm.
#---------------------------------------------------------------------------
def get_extract_files(data_dir, model, sdt_yyyymm, edt_yyyymm):
    files_out = []
    files = glob.glob(data_dir + '/' + model + '*' + '.nc')
    for f in range(len(files)):
        dt_c = re.findall(r'\.(\d{10})\.', files[f])
        if (dt_c):
            dt_c = dt_c[0]
            yyyymm_c = dt_c[0:6]
            if (yyyymm_c >= sdt_yyyymm and yyyymm_c <= edt_yyyymm):
                files_out.append(files[f])
    return sorted(files_out)

#---------------------------------------------------------------------------
# Get seasonal stats.
#---------------------------------------------------------------------------
def get_seasonal_files(files):

    springs = defaultdict(list)
    summers = defaultdict(list)
    falls   = defaultdict(list)
    winters = defaultdict(list)
    spring_ndays = {}
    summer_ndays = {}
    fall_ndays   = {}
    winter_ndays = {}
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

    return springs, summers, falls, winters, years

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
    
    for s in range(len(stns)):
        stn = stns[s]
        for m in range(len(models)):
            mod = models[m]
            ndays = {}
            years_mod_stn = []
            for d in range(len(dts_all)):
                yyyy = dts_all[d][0:4]
                mm = dts_all[d][4:6]
                years_mod_stn.append(yyyy)

                if ((mod,var,stn,dts_all[d]) in data_all):
                    dat_c = data_all[mod,var,stn,dts_all[d]]
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
    plt.ylabel(ylabs[var+stat], fontsize=fs+1)
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
    xlabs = []
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
                        ptot.append(statplot[key] * mm2in)
                    elif (var == 'T2MAX' or var == 'T2MIN'):
                        ptot.append(statplot[key] - 273.15)
                    elif (var.find('SPD') >= 0):
                        ptot.append(statplot[key])
                    if (m == 0):
                        xlabs.append(season.title() + ' ' + years[y])
        plt.plot(ptot, alpha = 0.3, label = mod.upper(), color = cols[m])
        lab_c = mod.upper() + ', ' + str(smooth_fact) + '-pt smoothing'
        plt.plot(smooth(ptot, smooth_fact), label = lab_c, color = cols[m])

    otot = []
    for y in range(len(years_obs)):
        yyyy = years[y]
        for ss in range(len(seasons)):
            season = seasons[ss]
            key = (season,stn,yyyy)
            if (season in seasonsplot):
                otot.append(obsplot[key])
    plt.plot(otot, alpha=0.3, label='Observed', color='black')
    lab_c = 'Observed ' + str(smooth_fact) + '-pt smoothing'
    plt.plot(smooth(otot,smooth_fact), label=lab_c, color='black')

    #--- y-axis labeling.
    plt.ylim(ylims[var+stat])
    plt.ylabel(ylabs[var+stat], fontsize=fs+1)
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
def mapper(var, lat, lon, grid, levs, cmap, maplw, colorbar_lab, titlein, \
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

#    print 'cmap = ', cmap
#    print 'cmaps[var] = ', cmaps[var]
#    print 'cmap_temp = ', cmap_temp
#    cs = plt.contourf(x, y, grid, levs, cmap=cmaps[var])
    cs = plt.contourf(x, y, grid, levs, cmap=cmap)    

    cbar = map.colorbar(cs, location='bottom', pad="3%", size=0.1, ticks=levs)
    cbar.set_label(colorbar_lab, fontsize=fs, size=fs-1)
    cbar.ax.tick_params(labelsize=fs-1)

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
