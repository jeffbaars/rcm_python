#!/usr/bin/python
import os, sys, glob, re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid

#from netCDF4 import Dataset as NetCDFFile
#from collections import defaultdict
#from utils_cmap import *
#import pickle
#from utils_date import *



def plot_it(data_all, dts_all, models, stns, var, stat, cols, \
            data_obs, dts_obs):

    mod = 'gfdl-cm3'
    var = 'T2MAX'
    stn = 'KSEA'
    stat = 'avg'
    plotfname = 'test.png'

    fig, ax = plt.subplots( figsize=(10,6) )
    lw_sm = 1.8

    for m in range(len(models)):
        ptot = []
        mod = models[m]
        for d in range(len(dts_all)):
            dt_c = dts_all[d]
            yyyy = dt_c[0:4]
            mm = dt_c[4:6]
            if (yyyy == '1999' and (mm >= '06' and mm <= '08')):
                key = (mod,var,stn,dt_c)
                if (key in data_all):
                    mod_c = data_all[key]
                    ptot.append(mod_c - 273.15)
        lab_c = mod.upper()
        plt.plot(ptot, label = lab_c, color = cols[m], linewidth=lw_sm)

    ptot = []
    for d in range(len(dts_obs)):
        dat_c = data_obs['TMAX','KSEA',dts_obs[d]]
        ptot.append(dat_c)
    lab_c = 'Obs'
    plt.plot(ptot, label = lab_c, color = 'black', linewidth=2.1)

    print 'xli ', plotfname, ' &'
    plt.savefig(plotfname)
    plt.close()

    return -9999
