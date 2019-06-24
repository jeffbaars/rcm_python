#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
#from make_cmap import *
#from netCDF4 import Dataset as NetCDFFile
#from utils_wrfout import *
#from utils_date import *
#import pickle

#from ghcndextractor import *
#from ghcndextractor import ghcndextractor
import codecs

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir    = '/home/disk/spock/jbaars/rcm'
py_dir     = rcm_dir + '/python'
plot_dir   = rcm_dir + '/plots'
pickle_dir = rcm_dir + '/pickle'
#ghcnFolder  = rcm_dir + '/obs'
ghcnd_dir   = rcm_dir + '/obs/ghcnd_all'
data_dir   = '/home/disk/r2d2/steed/cmip5/rcp8.5'

geo_em = '/home/disk/a125/steed/run/geo_em.d02.nc'

#---------------------------------------------------------------------------
# Settings.
#---------------------------------------------------------------------------
vars = ['PREC', 'T2MAX', 'T2MIN', 'SPDUV10MAX', 'SPDUV10MEAN', 'SPDUV10STD']
#sdt = '202001'
#edt = '202912'
sdt = '197001'
edt = '209912'

dtfmt = "%Y%m%d%H"

stns   = ['KSEA', 'KSMP', 'KYKM', 'KGEG', 'KPDX', 'KMFR', 'KHQM']
models = ['gfdl-cm3', 'miroc5']
stn_cols = ['b', 'r', 'g', 'c', 'm', 'k', 'grey']
mod_cols = ['b', 'r']
latpts = [47.44472, 47.27667, 46.56417, 47.62139, 45.59083, \
          42.38111, 46.97278]
lonpts = [-122.31361, -121.33722, -120.53361, -117.52778, -122.60028, \
          -122.87222, -123.93028]
elevs  = [130.0, 1207.0, 333.0, 735.0, 8.0, 405.0, 4.0]

#---------------------------------------------------------------------------
#
#---------------------------------------------------------------------------
stn_ghcnd = 'USW00024233'
file_c = ghcnd_dir + '/' + stn_ghcnd + '.dly'

if not os.path.isfile(file_c):
    print 'i do not see ', file_c
    sys.exit()

print 'i see ', file_c
readLoc = codecs.open(file_c, "r", "utf-8")
allLines = readLoc.readlines()
readLoc.close
    
for lineOfData in allLines:
    print lineOfData
    countryCode = lineOfData[0:2]
    stationID = lineOfData[0:11]
    stationMonthCode = lineOfData[0:17]
    year = lineOfData[11:15]
    month = lineOfData[15:17]
    element = lineOfData[17:21]

#    print 'year, month, element = ', year, month, element
    for x in range(0, 30):
        dayOM = x + 1
        offsetStart = (x*8)+21
        offsetEnd = offsetStart + 8
        ld = lineOfData[offsetStart:offsetEnd]
        val = ld[0:5]
        mflag = ld[5:6]
        qflag = ld[6:7]
        sflag = ld[7:8]

#        print '-----'
#        print 'ld = ', element, ld
#        print 'val = ', val
#        print 'mflag = ', mflag
#        print 'qflag = ', qflag
#        print 'sflag = ', sflag        

        if val == -9999.0:
            val = np.nan
        
        if element == 'TMAX' or element == 'TMIN':
            #--- Temperatures are in tenths of degrees C.
            valout = float(val) / 10.0
        elif element == 'PRCP':
            #--- Precip is in tenths of a mm.  Converting to inches.
            valout = (float(val) / 10.0) * 0.0393701
            
        if qflag and qflag.strip():
            print 'qflag is not empty, mflag, qflag, sflag = ', \
                  mflag, ',', qflag, ',', sflag
            valout = np.nan            

        print element, year, month, valout
#        sys.exit()
    
#    sys.exit()
        


#readLoc = codecs.open(file_c, "r", "utf-8")
#allLines = readLoc.readlines()
#readLoc.close
#
#for eachReadLine in allLines:
#    readRow(eachReadLine)


#ghcndextractor.countries = ["US"]
#ghcndextractor.states = ["NJ"]     
#ghcndextractor.ghcnFolder = "/Users/d035331/Documents/Demo/DemoData_Weather/ghcnd_all"   
#     
#ghcndextractor.getStationsFromFiles()
#stations = ghcndextractor.getCSVStationMetaData()
#ghcndextractor.readDailyFiles()
#dayCSV = ghcndextractor.getDailyDataCSV(["12"], ["25"], ["USW00014780"])
#print(dayCSV)
#
#
#
#
#
#stn_ghcnd = 'USW00024233'
#
##file_c = ghcnd_dir + '/' + stn_ghcnd + '.dly'
#
#stationIDCodes = ['USW00024233']
#iret = readDailyFiles()


sys.exit()

