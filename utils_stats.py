#!/usr/bin/python
import os, sys, glob, re, math
import numpy as np
from collections import defaultdict
import pickle
from utils_date import *
from utils_plot import *
from utils_ghcnd_obs import *

dtfmt = "%Y%m%d%H"

mm2in           = 0.0393700787
min_season_days = 85
std_lapse       = 0.0065  #-- std. atmos. lapse rate
mvc             = -9999.0

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
def get_ensmean_std(var, models, years, season, \
                    sum_files, spr_files, win_files, fal_files, ann_files, \
                    pickle_file):
    nyears = 0
    nmodels = len(models)
    files = []
    dts_all = []
    for y in range(len(years)):
        yyyy = years[y]

        for m in range(nmodels):
            model = models[m]
            sum_files_c = sum_files[model]
            spr_files_c = spr_files[model]
            win_files_c = win_files[model]
            fal_files_c = fal_files[model]
            ann_files_c = ann_files[model]

            if season == 'summer':
                if (sum_files_c[yyyy] == sum_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(sum_files_c[yyyy])):
                        files.append(sum_files_c[yyyy][f])
                    else:
                        continue
            elif season == 'spring':
                if (spr_files_c[yyyy] == spr_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(spr_files_c[yyyy])):
                        files.append(spr_files_c[yyyy][f])
                    else:
                        continue
            elif season == 'winter':
                if (win_files_c[yyyy] == win_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(win_files_c[yyyy])):
                        files.append(win_files_c[yyyy][f])
                    else:
                        continue
            elif season == 'fall':                        
                if (fal_files_c[yyyy] == fal_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(fal_files_c[yyyy])):
                        files.append(fal_files_c[yyyy][f])
                    else:
                        continue
            elif season == 'annual':
                if (ann_files_c[yyyy] == ann_files_c[yyyy]):
                    nyears += 1
                    for f in range(len(ann_files_c[yyyy])):
                        files.append(ann_files_c[yyyy][f])
                    else:
                        continue

    print 'files = ', files
    
    #--- Get sorted, unique set of dates from files found above.
    nfiles = len(files)
    dts_all = []
    for f in range(nfiles):
        dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]        
        dts_all.append(dt_c)
    dts_all = list(sorted(set(dts_all)))

    #--- For each unique date, read in files for each model.
    for d in range(len(dts_all)):

        #--- First make sure all members have a file.
        files_c = []
        print dts_all[d]
        for f in range(nfiles):
            dt_c = re.findall(r'\.(\d{10})\.', files[f])[0]
            if dt_c == dts_all[d]:
                files_c.append(files[f])
        if len(files_c) != nmodels:
            print 'missing member(s) for date = ', dts_all[d]
            sys.exit()

        #--- Now read in files for each member.
        for n in range(len(files_c)):
            print 'loading ', files_c[n]
            (vardat, ntimes, dt) = load_wrfout(files_c[n], [var])
            plotdat = vardat[var][:,:,:]
            nx,ny,nt = plotdat.shape
            if n == 0:
                dat_all = np.ones((nx, ny, nt, nmodels)) * np.nan
                stdevs_tmp = np.ones((nx, ny, nt)) * np.nan                
            if n == 0 and d == 0:
                stdevs = np.ones((nx, ny, nt)) * np.nan
            dat_all[:,:,:,n] = plotdat

        #--- For each ntimes (these should be monthly files, so ntimes are
        #--- days and will be ~30 days), calculate a standard deviation.
        for nt in range(ntimes):
            dat_c = np.squeeze(dat_all[:,:,nt,:])
            stdevs_tmp[:,:,nt] = np.std(dat_c, axis=2)

        #--- Add this month's standard deviation grids to our final array,
        #--- stdevs.
        if d == 0:
            stdevs = stdevs_tmp
        else:
            stdevs = np.concatenate((stdevs, stdevs_tmp), axis = 2)

    #--- Calculate the mean of the daily standard deviations.
    stdevmean = np.mean(stdevs, axis=2)
        
    #--- Create pickle file.
    pickle.dump((stdevmean), open(pickle_file,'wb'), -1)

    return stdevmean

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
# Get seasonal stats (totals or averages currently).
#---------------------------------------------------------------------------
def get_seasonal_stats(data_all, dts_all, models, stns, var, stat):

    out_sum   = {}
    out_avg   = {}
    out_max   = {}
    out_min   = {}
    years_all = []
    
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

                #--- Very first date of all model runs will have unusable
                #--- (all 0's) in the fields.  So don't try and use it.
                if dts_all[d] == '1970010100':
                    continue

                if ((var,stn,dts_all[d]) in data_all_c):
                    dat_c = data_all_c[var,stn,dts_all[d]]
                else:
                    ndays[season+yyyy] = np.nan
                    out_sum[season,mod,stn,yyyy] = np.nan
                    out_max[season,mod,stn,yyyy] = np.nan
                    out_min[season,mod,stn,yyyy] = np.nan
                    continue

                #--- annual is all days so get sum/max/min for every dts_all.
                (ndays, out_sum, out_max, out_min) = \
                        get_season_summaxmin(dat_c, mod, stn, 'annual',yyyy,mm,\
                                             ndays, out_sum, out_max, out_min)

                if (mm == '03' or mm == '04' or mm == '05'):
                    season = 'spring'
                if (mm == '06' or mm == '07' or mm == '08'):
                    season = 'summer'
                if (mm == '09' or mm == '10' or mm == '11'):
                    season = 'fall'
                if (mm == '12' or mm == '01' or mm == '02'):
                    season = 'winter'
                (ndays, out_sum, out_max, out_min) = \
                        get_season_summaxmin(dat_c, mod, stn, season, yyyy, mm,\
                                             ndays, out_sum, out_max, out_min)
                        
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

#----------------------------------------------------------------------------
# Add data to sums, max's and min's arrays.
#----------------------------------------------------------------------------
def get_season_summaxmin(dat_c, mod, stn, season, yyyyin, month, ndays, \
                         out_sum, out_max, out_min):

    #--- Look out for cases where data is NaN (e.g. the 31st of June) and
    #--- do nothing (return) in those cases.
    if np.isnan(dat_c):
        return ndays, out_sum, out_max, out_min
        
    #--- Put December winter data into next year's winter-- naming
    #--- winters by their Jan/Feb year rather than their Dec year.
    if month == '12':
        yyyy = str(int(yyyyin) + 1)
    else:
        yyyy = yyyyin
        
    #--- Count number of days for this season and year.
    if (season+yyyy) in ndays:
        ndays[season+yyyy] += 1
    else:
        ndays[season+yyyy] = 1

    key = (season,mod,stn,yyyy)

    #--- Get sum of sent in data.
    if key in out_sum:
        out_sum[key] = out_sum[key] + dat_c
    else:
        out_sum[key] = dat_c

    #--- Get maximum.
    if key in out_max:
        if dat_c > out_max[key]:
            out_max[key] = dat_c
    else:
        out_max[key] = dat_c
            
    #--- Get minimum.
    if key in out_min:
        if dat_c < out_min[key]:
            out_min[key] = dat_c
    else:
        out_min[key] = dat_c

    return ndays, out_sum, out_max, out_min
