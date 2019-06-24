#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt

#from datetime import date, timedelta, datetime
#import time, stat
#from collections import defaultdict
#from utils_date import *
#from mpl_toolkits.basemap import Basemap, cm


mm2in = 0.0393700787
min_season_days = 85

fs      = 9
titlefs = 9
width   = 10
height  = 6

#---------------------------------------------------------------------------
# Load geo_em file.
#---------------------------------------------------------------------------
def load_geo_em(geo_em):
    nc = NetCDFFile(geo_em, 'r')
    lat = nc.variables['XLAT_M'][0,:,:]
    lon = nc.variables['XLONG_M'][0,:,:]
    lat = np.transpose(lat)
    lon = np.transpose(lon)
    nc.close()
    return lat, lon

#---------------------------------------------------------------------------
# Load wrfout file.
#---------------------------------------------------------------------------
def load_wrfout(wrfout, var):
    nc = NetCDFFile(wrfout, 'r')
    ntimes = len(nc.dimensions['Time'])
    dataout = nc.variables[var][:,:,:]
    dataout = np.transpose(dataout)

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

    return files_out

#---------------------------------------------------------------------------
# Get seasonal totals.
#---------------------------------------------------------------------------
def get_seasonal_totals(data_all, dts_all, models, stns):

    springs = {}
    summers = {}
    falls   = {}
    winters = {}
    years_all = []
    
    for s in range(len(stns)):
        stn = stns[s]
        for m in range(len(models)):
            mod = models[m]

            spring_ndays = {}
            summer_ndays = {}
            fall_ndays   = {}
            winter_ndays = {}
            years_mod_stn = []
            
            for d in range(len(dts_all)):
                yyyy = dts_all[d][0:4]
                mm = dts_all[d][4:6]
                years_mod_stn.append(yyyy)

                #--- Spring.
                if (mm == '03' or mm == '04' or mm == '05'):

                    #--- Count number of days for this year's spring.
                    if ((yyyy) in spring_ndays):
                        spring_ndays[yyyy] += 1
                    else:
                        spring_ndays[yyyy] = 1

                    if ((mod,stn,yyyy) in springs):
                        springs[mod,stn,yyyy] = springs[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        springs[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]
                        
                #--- Summer.
                if (mm == '06' or mm == '07' or mm == '08'):
                    #--- Count number of days for this year's summer.
                    if ((yyyy) in summer_ndays):
                        summer_ndays[yyyy] += 1
                    else:
                        summer_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in summers):
                        summers[mod,stn,yyyy] = summers[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        summers[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]

                #--- Fall.
                if (mm == '09' or mm == '10' or mm == '11'):
                    #--- Count number of days for this year's fall.
                    if ((yyyy) in fall_ndays):
                        fall_ndays[yyyy] += 1
                    else:
                        fall_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in falls):
                        falls[mod,stn,yyyy] = falls[mod,stn,yyyy] + \
                                              data_all[mod,stn,dts_all[d]]
                    else:
                        falls[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]

                #--- Winter, December only.  Put December winter data into
                #--- next year's winter-- naming winters by their Jan/Feb
                #--- year rather than their Dec year.
                if (mm == '12'):
                    yyyyw = str(int(yyyy) + 1)
                    #--- Count number of days for this year's winter.
                    if ((yyyyw) in winter_ndays):
                        winter_ndays[yyyyw] += 1
                    else:
                        winter_ndays[yyyyw] = 1
                    if ((mod,stn,yyyyw) in winters):
                        winters[mod,stn,yyyyw] = winters[mod,stn,yyyyw] + \
                                                 data_all[mod,stn,dts_all[d]]
                    else:
                        winters[mod,stn,yyyyw] = data_all[mod,stn,dts_all[d]]
                #--- Winter, Jan & Feb.
                if (mm == '01' or mm == '02'):
                    #--- Count number of days for this year's winter.
                    if ((yyyy) in winter_ndays):
                        winter_ndays[yyyy] += 1
                    else:
                        winter_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in winters):
                        winters[mod,stn,yyyy] = winters[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        winters[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]

            years = list(sorted(set(years_mod_stn)))

            for y in range(len(years)):
                yyyy = years[y]
                if (spring_ndays[yyyy] < min_season_days):
                    springs[mod,stn,yyyy] = np.nan
                if (summer_ndays[yyyy] < min_season_days):
                    summers[mod,stn,yyyy] = np.nan
                if (fall_ndays[yyyy] < min_season_days):
                    falls[mod,stn,yyyy] = np.nan
                if (winter_ndays[yyyy] < min_season_days):
                    winters[mod,stn,yyyy] = np.nan

            years_all.append(years_mod_stn)

    years = list(sorted(set(years_all[0])))            
            
    return springs, summers, falls, winters, years

#---------------------------------------------------------------------------
# Make time series for seasonsal total precip for a set of stations.
#---------------------------------------------------------------------------
def ts(winters, springs, summers, falls, years, stns, mod, titlein, plotfname):

    fig, ax = plt.subplots( figsize=(10,6) )
    title1 = mod + ', ' + str(years[0]) + ' - ' + str(years[len(years)-1])

    xlabs = []
    for s in range(len(stns)):
        ptot = []
        stn = stns[s]
        print stn
        for y in range(len(years)):
            yyyy = years[y]
            ptot.append(winters[mod,stns[s],yyyy] * mm2in)
            ptot.append(springs[mod,stns[s],yyyy] * mm2in)
            ptot.append(summers[mod,stns[s],yyyy] * mm2in)
            ptot.append(falls[mod,stns[s],yyyy]   * mm2in)
            if (s == 0):
                xlabs.append('Winter ' + years[y])
                xlabs.append('Spring ' + years[y])
                xlabs.append('Summer ' + years[y])
                xlabs.append('Fall '   + years[y])            
        print stn, ptot
        plt.plot(ptot, label = stn)

    plt.ylim([0, 40.0])
    plt.xlabel('Season', fontsize=fs+1)
    plt.ylabel('Forecast Seasonal Total Precipitation (in.)', fontsize=fs+1)

    plt.xticks(range(0,len(xlabs)))
    plt.tick_params(axis='x', which='major', labelsize=fs-1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    
    plt.title(title1, fontsize=titlefs, fontweight='bold')

    ax.set_xticklabels( xlabs, rotation=90)
    plt.tight_layout()
    plt.grid()
    plt.legend(fontsize = fs, loc = 'best')

    print 'creating ', plotfname
    plt.savefig(plotfname)
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
