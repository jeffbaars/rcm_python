#!/usr/bin/python
import os, sys

levs_temp = [ 254, 256, 258, 260, 262, 264, 266, 268, 270, 272, 274, 276, \
              278, 280, 282, 284, 286, 288, 290, 292, 294, 296, 298, 300, \
              302, 304, 306, 308, 310, 312]
cmap_temp = [ (109, 227, 255), \
              (175, 240, 255), \
              (255, 196, 226), \
              (255, 153, 204), \
              (255,   0, 255), \
              (128,   0, 128), \
              (  0,   0, 128), \
              ( 70,  70, 255), \
              ( 51, 102, 255), \
              (133, 162, 255), \
              (255, 255, 255), \
              (204, 204, 204), \
              (179, 179, 179), \
              (153, 153, 153), \
              ( 96,  96,  96), \
              (128, 128,   0), \
              (  0,  92,   0), \
              (  0, 128,   0), \
              ( 51, 153, 102), \
              (157, 213,   0), \
              (212, 255,  91), \
              (255, 255,   0), \
              (255, 184, 112), \
              (255, 153,   0), \
              (255, 102,   0), \
              (255,   0,   0), \
              (188,  75,   0), \
              (171,   0,  56), \
              (128,   0,   0), \
              (163, 112, 255)]

levs_pcp = [ 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 
             27.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 80.0]
cmap_pcp = [(255,255,255), \
            (192,192,192), \
            (128,128,128 ), \
            (0,255,255), \
            (32,178,170), \
            (0,255,0), \
            (0,128,0), \
            (255,0,204), \
            (199,21,133), \
            (0,0,255), \
            (0,0,128), \
            (255,255,0), \
            (255,204,17), \
            (255,69,0), \
            (0,0,0),\
            (255,255,255)]
levs = { 'PREC':  levs_pcp,
         'T2MAX': levs_temp,
         'T2MIN': levs_temp
         }
cmaps = { 'PREC':  cmap_pcp,
          'T2MAX': cmap_temp,
          'T2MIN': cmap_temp
          }
