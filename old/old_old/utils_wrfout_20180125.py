#!/usr/bin/python
import os, sys, glob, re
import numpy as np
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from collections import defaultdict
from utils_cmap import *

seasons = ['spring', 'summer', 'fall', 'winter']

std_lapse = 0.0065  #-- std. atmos. lapse rate

mm2in = 0.0393700787
min_season_days = 85

fs      = 9
titlefs = 9
width   = 10
height  = 8
maplw   = 1.0

minlat =   39.0
maxlat =   51.0
minlon = -131.0
maxlon = -108.0

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
        if (spring_ndays[yyyy] < 3):
            springs[yyyy] = np.nan
        if (summer_ndays[yyyy] < 3):
            summers[yyyy] = np.nan
        if (fall_ndays[yyyy] < 3):
            falls[yyyy] = np.nan
        if (winter_ndays[yyyy] < 3):
            winters[yyyy] = np.nan

    return springs, summers, falls, winters, years

#---------------------------------------------------------------------------
# Get seasonal stats (totals or averages currently).
#---------------------------------------------------------------------------
def get_seasonal_stats(data_all, dts_all, models, stns, var, stat):

    springs = {}
    summers = {}
    falls   = {}
    winters = {}
    if (stat == 'avg'):
        springs_avg = {}
        summers_avg = {}
        falls_avg   = {}
        winters_avg = {}
    years_all = []

    stn_pr = 'KSEA'
    
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
                                                data_all[mod,var,stn,dts_all[d]]
                    else:
                        springs[mod,stn,yyyy] = data_all[mod,var,stn,dts_all[d]]

#                    if (stn == stn_pr):
#                        print '-----'
#                        print s, stn, m, mod, spring_ndays[yyyy], yyyy
#                        print dts_all[d], data_all[mod,var,stn,dts_all[d]]
#                        print 'springs[mod,stn,yyyy] = ', springs[mod,stn,yyyy]
                        
                #--- Summer.
                if (mm == '06' or mm == '07' or mm == '08'):
                    #--- Count number of days for this year's summer.
                    if ((yyyy) in summer_ndays):
                        summer_ndays[yyyy] += 1
                    else:
                        summer_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in summers):
                        summers[mod,stn,yyyy] = summers[mod,stn,yyyy] + \
                                                data_all[mod,var,stn,dts_all[d]]
                    else:
                        summers[mod,stn,yyyy] = data_all[mod,var,stn,dts_all[d]]

                #--- Fall.
                if (mm == '09' or mm == '10' or mm == '11'):
                    #--- Count number of days for this year's fall.
                    if ((yyyy) in fall_ndays):
                        fall_ndays[yyyy] += 1
                    else:
                        fall_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in falls):
                        falls[mod,stn,yyyy] = falls[mod,stn,yyyy] + \
                                              data_all[mod,var,stn,dts_all[d]]
                    else:
                        falls[mod,stn,yyyy] = data_all[mod,var,stn,dts_all[d]]

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
                                                 data_all[mod,var,stn,\
                                                          dts_all[d]]
                    else:
                        winters[mod,stn,yyyyw] = data_all[mod,var,stn,\
                                                          dts_all[d]]
                #--- Winter, Jan & Feb.
                if (mm == '01' or mm == '02'):
                    #--- Count number of days for this year's winter.
                    if ((yyyy) in winter_ndays):
                        winter_ndays[yyyy] += 1
                    else:
                        winter_ndays[yyyy] = 1
                    if ((mod,stn,yyyy) in winters):
                        winters[mod,stn,yyyy] = winters[mod,stn,yyyy] + \
                                                data_all[mod,var,stn,dts_all[d]]
                    else:
                        winters[mod,stn,yyyy] = data_all[mod,var,stn,dts_all[d]]

            #--- Get list of unique, sorted years found above.
            years = list(sorted(set(years_mod_stn)))

            #--- Check number of data points per unique years and nan out
            #--- ones that do not have enough (< min_season_days).
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

            #--- if stat requested was 'avg', do averaging.
            if (stat == 'avg'):
                for y in range(len(years)):
                    yyyy = years[y]
                    springs_avg[mod,stn,yyyy] = springs[mod,stn,yyyy] / \
                                                spring_ndays[yyyy]

#                    print '-----'
#                    print stn,yyyy,springs[mod,stn,yyyy], spring_ndays[yyyy]
#                    print springs_avg[mod,stn,yyyy]
                    
                    summers_avg[mod,stn,yyyy] = summers[mod,stn,yyyy] / \
                                                summer_ndays[yyyy]
                    falls_avg[mod,stn,yyyy]   = falls[mod,stn,yyyy]   / \
                                                fall_ndays[yyyy]
                    winters_avg[mod,stn,yyyy] = winters[mod,stn,yyyy] / \
                                                winter_ndays[yyyy]
                
            #--- Keep around all years found for this model and station.
            years_all.append(years_mod_stn)

    #--- Final, unique, sorted list of years found.
    years = list(sorted(set(years_all[0])))            

    if (stat == 'avg'):
        return springs_avg, summers_avg, falls_avg, winters_avg, years
    else:
        return springs, summers, falls, winters, years

#---------------------------------------------------------------------------
# Make time series for ??? for a set of stations.
#---------------------------------------------------------------------------
def ts(winters, springs, summers, falls, years, stns, mod, \
       titlein, plotfname, var, cols):

    fig, ax = plt.subplots( figsize=(10,6) )
    xlabs = []
    for s in range(len(stns)):
        ptot = []
        stn = stns[s]
        for y in range(len(years)):
            yyyy = years[y]
            if (var == 'PREC'):
                ptot.append(winters[mod,stns[s],yyyy] * mm2in)
                ptot.append(springs[mod,stns[s],yyyy] * mm2in)
                ptot.append(summers[mod,stns[s],yyyy] * mm2in)
                ptot.append(falls[mod,stns[s],yyyy]   * mm2in)
            elif (var == 'T2MAX' or var == 'T2MIN'):
                ptot.append(winters[mod,stns[s],yyyy] - 273.15)
                ptot.append(springs[mod,stns[s],yyyy] - 273.15)
                ptot.append(summers[mod,stns[s],yyyy] - 273.15)
                ptot.append(falls[mod,stns[s],yyyy]   - 273.15)
            if (s == 0):
                xlabs.append('Winter ' + years[y])
                xlabs.append('Spring ' + years[y])
                xlabs.append('Summer ' + years[y])
                xlabs.append('Fall '   + years[y])            

#        window = 10
#        weights = np.repeat(1.0, window)/window
#        smas = np.convolve(values, weights, 'valid')

        plt.plot(ptot, alpha = 0.3, color = cols[s])
        plt.plot(smoothListGaussian(ptot, 3), label = stn, color = cols[s])

    #--- y-axis labeling.
    if (var == 'PREC'):
        ylab = 'Forecast Seasonal Total Precipitation (in.)'
        ylim1 = [0, 40.0]
    elif (var == 'T2MAX'):
        ylab = 'Forecast Seasonal Average Maximum Temperature ($^\circ$C)'
        ylim1 = [-10, 40.0]
    elif (var == 'T2MIN'):
        ylab = 'Forecast Seasonal Average Minimum Temperature ($^\circ$C)'
        ylim1 = [-10, 30.0]
    plt.ylim(ylim1)
    plt.ylabel(ylab, fontsize=fs+1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    

    #--- x-axis labels.
#    xticks_c = range(0,len(xlabs),10)
    xticks_c = np.arange(0,len(xlabs)-1,10)    
    print 'xticks_c = ', xticks_c
    
    plt.xticks(xticks_c)
    plt.tick_params(axis='x', which='major', labelsize=fs-1)
#    ax.set_xticklabels(xlabs, rotation=90)

#    print 'xlabs = ', xlabs
#    print 'xticks_c = ', xticks_c
#    print 'xlabs[0] = ', xlabs[[0,5,10]]
#    print 'xlabs[xticks_c] = ', xlabs[xticks_c]

#    ax.set_xticklabels(xlabs[xticks_c], rotation=90)
    ax.set_xticklabels(xlabs, rotation=90)        
    plt.xlabel('Season', fontsize=fs+1)

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    plt.tight_layout()
    #plt.grid()
    plt.legend(fontsize = fs, loc = 'best')

    print 'xli ', plotfname, ' &'
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

#---------------------------------------------------------------------------
# Mapper.
#---------------------------------------------------------------------------
def mapper(lat, lon, grid, minlat, maxlat, minlon, maxlon, \
           levs, cmap, maplw, colorbar_lab, titlein, plotfname):

    ur_lat = maxlat
    ur_lon = maxlon
    ll_lat = minlat
    ll_lon = minlon
    lat_ctr = ll_lat + ((ur_lat - ll_lat) * 0.5)
    lon_ctr = ll_lon + ((ur_lon - ll_lon) * 0.5)

    fig = plt.figure(figsize=(width,height))
    # left, bottom, width, height:
    ax = fig.add_axes([0.00,0.05,0.99,0.91])
    map = Basemap(resolution='i',projection='lcc',\
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

    cs = plt.contourf(x, y, grid, levs, cmap=cmap)

    cbar = map.colorbar(cs, location='bottom', pad="3%", size=0.1, ticks=levs)
    cbar.set_label(colorbar_lab, fontsize=fs, size=fs-1)
    cbar.ax.tick_params(labelsize=fs-1)

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    #--- Save plot.
    print 'creating ', plotfname
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
