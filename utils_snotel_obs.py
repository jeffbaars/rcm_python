#!/usr/bin/python
import os, sys, glob, re
import numpy as np
import csv

in2mm = 25.4

#---------------------------------------------------------------------------
# Read in a swe stations file.
#---------------------------------------------------------------------------
def read_swe_stations_file(station_file):
    stns      = []
    latpts    = []
    lonpts    = []
    elevs     = []
    stn_names = []
    with open(station_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            stns.append(row[0])
            latpts.append(float(row[1]))
            lonpts.append(float(row[2]))
            elevs.append(float(row[3]))
            s = row[4]
            s = s.strip(' "\'\t\r\n')
            stn_names.append(s)

    return stns, latpts, lonpts, elevs, stn_names

#---------------------------------------------------------------------------
# Read in a snotel data file. Columns are Water Year,Day,Oct,Nov,Dec,Jan,
# Feb,Mar,Apr,May,Jun,Jul,Aug,Sep.
#---------------------------------------------------------------------------
def read_snotel(swe_file):
    yyyy    = []
    swe_obs = []    
    with open(swe_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:

            #--- Skip comments.
            if row[0].startswith('#'):
                continue

            #--- Skip non-1st-of-the-month data.
            if row[1] != '01':
                continue

            #--- Skip if there's no actual value for this year.
            if not row[8]:
                continue

            #--- Grab April (1st).
            yyyy.append(row[0])
            swe_obs.append(float(row[8]) * in2mm)
    
    return yyyy, swe_obs

