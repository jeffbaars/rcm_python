#!/usr/bin/python
import os, sys

#--------------------------------------------------------------------------
# Rainwatch dBz and precip colors.
#--------------------------------------------------------------------------
# Color table taken from:
# http://pykl3radar.com/pykl3wiki/index.php/Color_table_customization
levs_dbz = [-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75]
cmap_dbz = [(255,255,255), \
            (255,255,255), \
            (4,233,227), \
            (4,158,243), \
            (4,2,243), \
            (2,250,2), \
            (17,121,1), \
            (0,140,0), \
            (255,255,0), \
            (228,187,2), \
            (255,148,0), \
            (252,0,0), \
            (210,0,0), \
            (187,0,0), \
            (247,0,252), \
            (151,84,197), \
            (255,255,255)]

#levs_pcp = [ 0.0, 0.01, 0.1, 0.2, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5,
#             3.0, 3.5, 4.0, 5.0, 6.0, 7.0 ]
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

levs_swe = [ 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
             1100, 1200, 1300, 1400, 1500, 1600]
cmap_swe = [(255,255,255), \
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

levs_swdown = [ 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, \
                600, 650, 700, 750, 800]
cmap_swdown = [(255,255,255), \
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

#--------------------------------------------------------------------------
# Temperature.
#--------------------------------------------------------------------------
#levs_temp = [ 270, 272, 274, 276, 278, 280, 282, 284, 286, \
#              288, 290]
#levs_temp = [ 264, 268, 272, 276, 280, 284, 288, 292, 296, 300, 304]
#cmap_temp = [( 5,   48,   97), \
#             (33,  102,  172), \
#             (67,  147,  195), \
#             (146, 197,  222), \
#             (209, 229,  240), \
#             (247, 247,  247), \
#             (254, 219,  199), \
#             (244, 165,  130), \
#             (214,  96,   77), \
#             (178,  24,   43), \
#             (103,   0,   31)]

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

#levs_temp_diff = [ -6, -4, -2, 0, 2, 4, 6, 8, 10, 12 ]
#levs_temp_diff = [ -3, -2, -1, 0, 1, 2, 3, 4, 6, 8 ]
#cmap_temp_diff = [ (  0,   0,  51), \
#                   (  0,   0, 204), \
#                   (153, 153, 255), \
#                   (255, 255, 255), \
#                   (255, 153, 153), \
#                   (255,  51,  51), \
#                   (204,   0,   0), \
#                   (102,   0,   0), \
#                   ( 51,   0,   0), \
#                   (  0,   0,   0)]
levs_temp_diff = [ -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8 ]
cmap_temp_diff = [ (  0,   0,  51), \
                   (  0,   0, 204), \
                   (153, 153, 255), \
                   (255, 255, 255), \
                   (255, 100, 100), \
                   (255, 153, 153), \
                   (255,  51,  51), \
                   (200,   0,   0), \
                   (150,   0,   0), \
                   (100,   0,   0), \
                   ( 50,   0,   0), \
                   (  0,   0,   0)]

levs_pcp = [ -600, -500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500, \
              600 ]
cmap_pcp = [ (128,   0, 128), \
             (  0,   0,  51), \
             (  0,   0, 204), \
             (153, 153, 255), \
             (255, 255, 255), \
             (255, 100, 100), \
             (255, 153, 153), \
             (255,  51,  51), \
             (200,   0,   0), \
             (150,   0,   0), \
             (100,   0,   0), \
             ( 50,   0,   0), \
             (  0,   0,   0)]

#levs_pcp = [ -500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500 ]
#levs_pcp = [ -1000, -800, -600, -400, -200, 0, 200, 400, 600, 800, 1000 ]
#levs_pcp = [ -1200, -800, -400, -200, -100, 0, 100, 200, 400, 800, 1200, 1600 ]
#levs_pcp = [ -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60 ]
levs_pcp = [ -50, -40, -20, -15, -10, -5, 0, 5, 10, 15, 20, 40, 50 ]
cmap_pcp = [( 63,  37, 11), \
            ( 84,  48,  5), \
            (140,  81, 10), \
            (191, 129, 45), \
            (223, 194,125), \
            (246, 232,195), \
            (245, 245,245), \
            (199, 234,229), \
            (128, 205,193), \
            ( 53, 151, 43), \
            (  1, 102, 95), \
            (  0,  60, 48), \
            (  0, 100,  0)]

#cmap_pcp = [( 84,  48,  5), \
#            (140,  81, 10), \
#            (191, 129, 45), \
#            (223, 194,125), \
#            (246, 232,195), \
#            (245, 245,245), \
#            (199, 234,229), \
#            (128, 205,193), \
#            ( 53, 151, 43), \
#            (  1, 102, 95), \
#            (  0,  60, 48)]


#levs_pcp_diff = [-500, -250, 0, 250, 500, 750]
#cmap_pcp_diff = [(  0,  0, 51), \
#                 (  0,  0,204), \
#                 (153,153,255), \
#                 (255,255,255), \
#                 (255,100,100), \
#                 (255,153,153)]

#levs_pcp_diff = [ 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, \
#             27.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 80.0]
#cmap_pcp_diff = [(255,255,255), \
#            (192,192,192), \
#            (128,128,128 ), \
#            (0,255,255), \
#            (32,178,170), \
#            (0,255,0), \
#            (0,128,0), \
#            (255,0,204), \
#            (199,21,133), \
#            (0,0,255), \
#            (0,0,128), \
#            (255,255,0), \
#            (255,204,17), \
#            (255,69,0), \
#            (0,0,0),\
#            (255,255,255)]
