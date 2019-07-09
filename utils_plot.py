#!/usr/bin/python
import os, sys, glob, re, math
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from utils_cmap import *
from utils_date import *

zooms = ['Z1', 'Z2']

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
    'T2MEAN': 'Average Mean Temperature ($^\circ$C)',    
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
    'T2MEANmin': [-30.0, 25.0],
    'T2MEANmax': [-10.0, 50.0],
    'T2MEANavg': [-10.0, 50.0],
    'SPDUV10MEANavg': [0, 5.0],
    'SPDUV10MAXmax': [0, 20.0],
    'SNOWtot': [0, 2000.0]
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
    'T2MEANavg': 'Average Mean Temperature ($^\circ$C)',
    'T2MEANmax': 'Max Mean Temperature ($^\circ$C)',
    'T2MEANmin': 'Min Mean Temperature ($^\circ$C)',
    'SPDUV10MEANavg': 'Average Wind Speed (m/s)',
    'SPDUV10MAXmax': 'Maximum Wind Speed (m/s)',
    'SNOWtot': 'Snow Water Equivalent (mm)'
    }

colorbar_labs = {
    'T2MAX': 'Temperature ($^\circ$C)',
    'T2MIN': 'Temperature ($^\circ$C)',
    'T2MEAN': 'Temperature ($^\circ$C)',    
    'PREC': 'Precipitation (in.)',
    'SPDUV10MEAN': 'Wind Speed (m/s)',
    'SPDUV10MAX': 'Wind Speed (m/s)',
    'SNOW': 'Snow Water Equivalent (mm)',
    'SWDOWN': 'Shortwave Down Radiation (W/m^2)'     
    }

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

fs      = 9
titlefs = 9
width   = 10
height  = 8
maplw   = 1.0

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
                    elif (var == 'T2MAX' or var == 'T2MIN' or var == 'T2MEAN'):
                        ptot.append(statplot[season,mod,stns[s],yyyy] - 273.15)
                    elif (var.find('SPD') >= 0):
                        ptot.append(statplot[season,mod,stns[s],yyyy])
                    if (s == 0):
                        xlabs.append(season.title() + ' ' + years[y])
        plt.plot(ptot, alpha = 0.3, label = stn, color = cols[s])
        plt.plot(ptot, label = stn, color = cols[s])        

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
                    elif (var == 'T2MAX' or var == 'T2MIN' or var == 'T2MEAN'):
                        ptot_c = statplot[key] - 273.15
                    elif (var.find('SPD') >= 0):
                        ptot_c = statplot[key]
                    if (m == 0):
                        xlabs.append(years[y])

                    ptot.append(ptot_c)
            mod_all[m,y] = ptot_c
            
        plt.plot(ptot, label = mod.upper(), \
                 color = cols[m], linewidth=lw_sm, alpha = 0.23) 

    #--- Plot ensemble mean.
    plt.plot(np.nanmean(mod_all, axis=0), label = 'Ensemble Mean', \
             color = 'darkgreen', linewidth=2.2)
    
    #--- Plot observations.
    otot = []
    for y in range(len(years_obs)):
        yyyy = years_obs[y]
        for ss in range(len(seasons)):
            season = seasons[ss]
            key = (season,stn,yyyy)
            if (season in seasonsplot):
                otot.append(obsplot[key])
    plt.plot(otot, label='Observed', color='black', \
             linestyle='None', marker='o', markersize=3)

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
# Make time series of sent-in data for a set of models for one station.
#---------------------------------------------------------------------------
def ts_swe(swe_all, years, swe_obs, years_obs, stn, s, models, titlein, \
           plotfname, var, stat, cols):

    (nm,ns,nf) = swe_all.shape
    fig, ax = plt.subplots( figsize=(10,6) )
    lw_sm = 1.8

    #--- Plot models.
    xlabs = []
    mod_all = np.ones((len(models),len(years))) * np.nan
    for m in range(len(models)):
        mod = models[m]

        swe_plot = []
        for f in range(nf):
            swe_plot.append(swe_all[m,s,f])

        plt.plot(swe_plot, label = mod.upper(), color = cols[m], \
                           linewidth=lw_sm, alpha = 0.23)

    #--- Plot ensemble mean.
    plt.plot(np.nanmean(swe_all[:,s,:], axis=0), label = 'Ensemble Mean', \
             color = 'darkgreen', linewidth=2.2)

    #--- Plot observations.
    obs_plot      = []
    obs_yyyy_plot = []    
    for y in range(len(years)):
        if years[y] in years_obs:
            i = years_obs.index(years[y])
            obs_plot.append(swe_obs[i])
        else:
            obs_plot.append(np.nan)            
    plt.plot(obs_plot, label='Observed', color='black', \
             linestyle='None', marker='o', markersize=3)

    #--- y-axis labeling.
    plt.ylim(ylims[var+stat])
    plt.ylabel(labels[var+stat], fontsize=fs+1)
    plt.tick_params(axis='y', which='major', labelsize=fs+1)    

    #--- x-axis labels.
    xticks_c = range(0,len(years), 4)
    xlabs_c = [ years[i] for i in xticks_c ]    
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
# Run mapper.
#---------------------------------------------------------------------------
def run_mapper(var, model, rundir, season, files, yyyy, \
               lat, lon, levs_c, my_cmap ):

    file_plot = rundir + '/' + model + '_' + season + '_' + yyyy + '_' + \
                var + '.nc'
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

    syscom = 'rm -f ' + file_plot
    os.system(syscom)
    
    titlein = model.upper() + ', ' + var_lab[var] + ', ' + \
              seasons_lab[season] + ' ' + yyyy

    plotdat = vardat[var][:,:,0]
    if (var == 'PREC'):
        plotdat = plotdat * mm2in
    else:
        plotdat = plotdat - 273.15

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

    cs = plt.contourf(x, y, grid, levs, cmap=cmap_in)

    cbar = map.colorbar(cs, location='bottom', pad="3%", size=0.1, ticks=levs)
    cbar.set_label(colorbar_lab, fontsize=fs, size=fs-1)
    cbar.ax.tick_params(labelsize=fs-2)

    if var != 'SNOW':
        csl = plt.contour(x, y, grid, levs, colors = 'black', linewidths=0.8)
        plt.clabel(csl, inline=1, fontsize=10, fmt='%2.1f', fontweights='bold')

    plt.title(titlein, fontsize=titlefs, fontweight='bold')

    #--- Save plot.
    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)

    plt.close()

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
