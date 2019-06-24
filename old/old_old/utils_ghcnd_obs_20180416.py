#!/usr/bin/python
import os, sys, glob, re
import numpy as np
import codecs

seasons = ['spring', 'summer', 'fall', 'winter']
min_season_days = 85

ghcnd_dict = {
    'KSEA': 'USW00024233',
    'KSMP': 'USW00024237',
    'KYKM': 'USW00024243',
    'KGEG': 'USW00024157',
    'KPDX': 'USW00024229',
    'KMFR': 'USW00024225',
    'KHQM': 'USW00094225'}

#---------------------------------------------------------------------------
# Read in a GHCND obs file.
#---------------------------------------------------------------------------
def read_ghcnd(stns, elements, ghcnd_dir, syyyymm, eyyyymm):

    data_out = {}
    dts_all = []

    for s in range(len(stns)):
        stn = stns[s]
        stn_ghcnd = ghcnd_dict[stn]
        file_c = ghcnd_dir + '/' + stn_ghcnd + '.dly'

        if not os.path.isfile(file_c):
            print 'cannot find file ', file_c
            continue
        else:
            print 'Reading ', stn, ' file ', file_c

        readLoc = codecs.open(file_c, "r", "utf-8")
        allLines = readLoc.readlines()
        readLoc.close
    
        for lineOfData in allLines:
            countryCode = lineOfData[0:2]
            stationID = lineOfData[0:11]
            stationMonthCode = lineOfData[0:17]
            year = lineOfData[11:15]
            month = lineOfData[15:17]
            
            #--- Skip lines not in 'elements'.
            element = lineOfData[17:21]
            if not element in elements:
                continue

            #--- Skip lines with dates falling outside sent-in date ranges.
            yyyymm = year + month
            if (yyyymm < syyyymm or yyyymm > eyyyymm):
                continue

            for x in range(0, 31):
                dayOM = x + 1
                offsetStart = (x*8)+21
                offsetEnd = offsetStart + 8
                ld = lineOfData[offsetStart:offsetEnd]
                val = ld[0:5]
                mflag = ld[5:6]
                qflag = ld[6:7]
                sflag = ld[7:8]
                
                if float(val) == -9999.0:
                    val = np.nan
        
                if element == 'TMAX' or element == 'TMIN':
                    #--- Temperatures are in tenths of degrees C.
                    valout = float(val) / 10.0
                elif element == 'PRCP':
                    #--- Precip is in tenths of a mm.  Converting to inches.
                    valout = (float(val) / 10.0) * 0.0393701
                else:
                    sys.exit('read_ghcnd:  element not recognized!')

                #--- Look at qflag and set valout to nan if it's flagged.
                #--- Some of these out of range values looked legit though!
                if qflag and qflag.strip():
                    valout = np.nan            

                day_c = '{:02d}'.format(x+1)
                dt_c = year + month + day_c
#                print 'dt_c = ', dt_c
                dts_all.append(dt_c)
                data_out[element,stn,dt_c] = valout

        #--- Get unique list of dates.
        dts_all = list(sorted(set(dts_all)))
        
    return data_out, dts_all

#---------------------------------------------------------------------------
# Get seasonal stats (totals or averages currently).
#---------------------------------------------------------------------------
def get_seasonal_stats_ghcnd(data_all, dts_all, stns, var, stat):

    out_sum = {}
    out_avg = {}
    out_max = {}
    out_min = {}
        
    years_all = []

    for s in range(len(stns)):
        stn = stns[s]
        ndays = {}
        years_stn = []
        for d in range(len(dts_all)):
            yyyy = dts_all[d][0:4]
            mm = dts_all[d][4:6]
            years_stn.append(yyyy)
                
            if ((var,stn,dts_all[d]) in data_all):
                dat_c = data_all[var,stn,dts_all[d]]
            else:
                continue
            
            #--- Spring.
            if (mm == '03' or mm == '04' or mm == '05'):
                season = 'spring'
                get_season_summaxmin(dat_c, stn, season, yyyy, mm, ndays, \
                                     out_sum, out_max, out_min)
            #--- Summer.
            if (mm == '06' or mm == '07' or mm == '08'):
                season = 'summer'
                get_season_summaxmin(dat_c, stn, season, yyyy, mm, ndays, \
                                     out_sum, out_max, out_min)
            #--- Fall.
            if (mm == '09' or mm == '10' or mm == '11'):
                season = 'fall'
                get_season_summaxmin(dat_c, stn, season, yyyy, mm, ndays, \
                                     out_sum, out_max, out_min)

            #--- Winter.  Putting December winter data into
            #--- next year's winter-- naming winters by their Jan/Feb
            #--- year rather than their Dec year.
            if (mm == '12' or mm == '01' or mm == '02'):
                season = 'winter'
                get_season_summaxmin(dat_c, stn, season, yyyy, mm, ndays, \
                                         out_sum, out_max, out_min)
                        
        #--- Get list of unique, sorted years found above.
        years = list(sorted(set(years_stn)))

        #--- Check number of data points per unique years and nan out
        #--- ones that do not have enough (< min_season_days).
        for y in range(len(years)):
            yyyy = years[y]
            for s in range(len(seasons)):
                season = seasons[s]
                if ((season+yyyy) in ndays):
                    if (ndays[season+yyyy] < min_season_days):
                        out_sum[season,stn,yyyy] = np.nan
                        out_max[season,stn,yyyy] = np.nan
                        out_min[season,stn,yyyy] = np.nan
                else:
                    ndays[season+yyyy] = np.nan
                    out_sum[season,stn,yyyy] = np.nan
                    out_max[season,stn,yyyy] = np.nan
                    out_min[season,stn,yyyy] = np.nan

        #--- if stat requested was 'avg', do averaging.
        if (stat == 'avg' or stat == 'max' or stat == 'min'):
            for y in range(len(years)):
                yyyy = years[y]
                for s in range(len(seasons)):
                    season = seasons[s]
                    out_avg[season,stn,yyyy] = out_sum[\
                        season,stn,yyyy] / ndays[season+yyyy]
            
        #--- Keep around all years found for this station.
        years_all.append(years_stn)

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

#----------------------------------------------------------------------------
# Add data to sums, max's and min's arrays, for any month other than December.
#----------------------------------------------------------------------------
def get_season_summaxmin(dat_c, stn, season, yyyyin, month, ndays, \
                         out_sum, out_max, out_min):

    #--- Putg December winter data into next year's winter-- naming
    #--- winters by their Jan/Feb year rather than their Dec year.
    if (month == '12'):
        yyyy = str(int(yyyyin) + 1)
    else:
        yyyy = yyyyin
        
    #--- Count number of days for this year's spring.
    if (season+yyyy in ndays):
        ndays[season+yyyy] += 1
    else:
        ndays[season+yyyy] = 1
        
    if ((season,stn,yyyy) in out_sum):
        out_sum[season,stn,yyyy] = out_sum[season,stn,yyyy] + dat_c
    else:
        out_sum[season,stn,yyyy] = dat_c

    #--- Get maximum.
    if ((season,stn,yyyy) in out_max):
        if (dat_c > out_max[season,stn,yyyy]):
            out_max[season,stn,yyyy] = dat_c
        else:
            out_max[season,stn,yyyy] = dat_c
            
    #--- Get minimum.
    if ((season,stn,yyyy) in out_min):
        if (dat_c < out_min[season,stn,yyyy]):
            out_min[season,stn,yyyy] = dat_c
        else:
            out_min[season,stn,yyyy] = dat_c

    return ndays, out_sum, out_max, out_min

