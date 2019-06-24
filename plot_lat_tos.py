#!/usr/bin/python
import  sys, os, os.path, time, glob, re
import numpy as np
from igtg_check_utils import *
import matplotlib.pyplot as plt
#import pickle
from make_cmap import *

#---------------------------------------------------------------------------
# Paths.
#---------------------------------------------------------------------------
rcm_dir  = '/home/disk/spock/jbaars/rcm'
py_dir   = rcm_dir + '/python'
plot_dir = rcm_dir + '/plots'

x = range(140,161)
dat = [ 285.5343, 285.7745, 286.0928, 286.3722, 286.6160, 286.6160, \
        286.8224, 286.3380, 286.5203, 286.7137, 286.8461, 286.9745, \
        287.1689, 287.3819, 287.6205, 287.6205, 287.4508, 287.9211, \
        288.3215, 288.6165, 288.8190 ]
#x = range(0,len(dat))

plt.plot(x, dat, 'b')
plt.plot(x, dat, 'bo')
plt.grid()

plt.title('j = 92, i = 140 - 160')
plt.savefig('lat_tos_test.png')
plt.close()

sys.exit()

