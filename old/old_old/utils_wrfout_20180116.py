#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile

#from datetime import date, timedelta, datetime
#import time, stat
#from collections import defaultdict
#from utils_date import *
#from mpl_toolkits.basemap import Basemap, cm
#import matplotlib.pyplot as plt

mm2in = 0.0393700787

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
# Get seasonal totals.  NEED TO SORT OUT HOW TO HANDLE WINTERS AT THE
# EDGES (FIRST YEAR AND LAST) WHICH WON"T HAVE A FULL 3 MONTHS.
#---------------------------------------------------------------------------
def get_seasonal_totals(data_all, dts_all, models, stns):
    springs = {}
    summers = {}
    falls   = {}
    winters = {}
    
    for s in range(len(stns)):
        stn = stns[s]
        for m in range(len(models)):
            mod = models[m]
            for d in range(len(dts_all)):
                yyyy = dts_all[d][0:4]
                mm = dts_all[d][4:6]
                if (mm == '03' or mm == '04' or mm == '05'):

#                   if (stn == stn_pr and yyyy == yyyy_pr):
#                       print dts_all[d], mm, data_all[mod,stn,dts_all[d]]

                    if ((mod,stn,yyyy) in springs):
                        springs[mod,stn,yyyy] = springs[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        springs[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]

#                    if (stn == stn_pr and yyyy == yyyy_pr):
#                        print 'springs tot = ', springs[mod,stn,yyyy]
                        
                if (mm == '06' or mm == '07' or mm == '08'):
                    if ((mod,stn,yyyy) in summers):
                        summers[mod,stn,yyyy] = summers[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        summers[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]
                if (mm == '09' or mm == '10' or mm == '11'):
                    if ((mod,stn,yyyy) in falls):
                        falls[mod,stn,yyyy] = falls[mod,stn,yyyy] + \
                                              data_all[mod,stn,dts_all[d]]
                    else:
                        falls[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]
                if (mm == '12' or mm == '01' or mm == '02'):
                    if ((mod,stn,yyyy) in winters):
                        winters[mod,stn,yyyy] = winters[mod,stn,yyyy] + \
                                                data_all[mod,stn,dts_all[d]]
                    else:
                        winters[mod,stn,yyyy] = data_all[mod,stn,dts_all[d]]

    return springs, summers, falls, winters


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
