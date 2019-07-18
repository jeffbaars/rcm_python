#!/usr/bin/python
#----------------------------------------------------------------------------
# Climate extreme indices based on Rick's code which is based on:
# http://etccdi.pacificclimate.org/list_27_indices.shtml
#----------------------------------------------------------------------------
import os, sys, glob, re, math
import numpy as np
from utils_date import *
from utils_ghcnd_obs import *

#basefirst = 1971
#baselast  = 2000
#first     = 1970
#last      = 2099

basefirst = 1971
baselast  = 1974
first     = 1970
last      = 1979

#--- percentages used for...
percs = [0.1, 0.9]

dpm = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

#----------------------------------------------------------------------------
# calculate daily percentages of the historical period.
#----------------------------------------------------------------------------
def get_daily_perc_mxmn(models, stns, mod_all_in, sdt):

    bdt = sdt + '0100'
    daily_perc = {}

    #--- Loop over models.
    for m in range(len(models)):
        model = models[m]
        mod_all = mod_all_in[model]
        #--- Loop over stations.
        for s in range(len(stns)):
            stn = stns[s]
            #--- Loop over each day of the year.
            for d in range(0,364):
                fdate = time_increment(bdt, d, dtfmt)
                basemmddhh = fdate[4:10]
                txunsorted = []
                tnunsorted = []
                #--- Loop over every year within the "base" period.
                for y in range(basefirst, baselast+1):
                    basedate = str(y) + basemmddhh
                    #--- Per CEI definition of daily percentage calc, Looks
                    #--- +/- 2 days from basedate, presumably to increase
                    #--- sample size(?).
                    for n in range(-2, 3):
                        dt = time_increment(basedate, n, dtfmt)
                        key = ('T2MAX', stn, dt)
                        if key in mod_all:
                            txunsorted.append(mod_all['T2MAX', stn, dt])
                        key = ('T2MIN', stn, dt)
                        if key in mod_all:
                            tnunsorted.append(mod_all['T2MIN', stn, dt])

                #--- Sort lists of maxes and mins.
                txsorted = list(sorted(txunsorted))
                tnsorted = list(sorted(tnunsorted))

                #--- Get percentage (percentages defined in percs) for
                #--- max and min.
                for p in range(len(percs)):
                    perc = percs[p]
                    np = int(round(float(len(txsorted)+1) * perc))
                    key = (model, 'T2MAX', stn, fdate, str(perc))
                    daily_perc[key] = txsorted[np]
                    
                    np = int(round(float(len(tnsorted)+1) * perc))
                    key = (model, 'T2MIN', stn, fdate, str(perc))
                    daily_perc[key] = tnsorted[np]
                    
    return daily_perc

#----------------------------------------------------------------------------
# Calculate monthly max and min climate index.
#----------------------------------------------------------------------------
def get_monthly_mxmn(daily_perc, mod_all, models, stns, sdt):

    model = models[0]
    stn = stns[0]

    for m in range(0,12):
        mstr = str(m+1).zfill(2)
        for y in range(first, last+1):
            bdt = str(y) + mstr + '0100'
            yyyymm = str(y) + mstr
            for d in range(0, dpm[m]):
                fdate = time_increment(bdt, d, dtfmt)
                basedate = '1970' + fdate[4:10]

                key   = ('T2MAX', stn, fdate)
                keydp = (model, 'T2MAX', stn, fdate)                

#                if mod_all[key] > daily_perc[keydp]

                print fdate
            sys.exit()


#    #--- Loop over models.
#    for m in range(len(models)):
#        model = models[m]
#        mod_all = mod_all_in[model]
#        #--- Loop over stations.
#        for s in range(len(stns)):
#            stn = stns[s]
#
#            for m in range(0,12):
#                print m
#                for y in range(first, last+1):
#                    print y

                    
        sys.exit()



