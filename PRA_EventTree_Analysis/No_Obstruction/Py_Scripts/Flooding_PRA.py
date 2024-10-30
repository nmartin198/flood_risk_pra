# -*- coding: utf-8 -*-
"""
.. module:: Flooding_PRA
   :platform: Windows, Linux
   :synopsis: No Obstruction branch, inundation simulations
 
.. moduleauthor:: Nick Martin <nick.martin@alumni.stanford.edu>

Provides a Monte Carlo implementation to execute all inundation simulations
for the No Obstruction branch of the event tree. This module is designed
for parallelization by running independent Monte Carlo simulations, by 
future weather realization, on different computers.

"""
# Copyright and License
"""
Copyright 2024 Vodanube LLC

Module Author: Nick Martin <nick.martin@alumni.stanford.edu>

This file is part of a Flood Risk PRA example study, hereafter Flood Risk PRA.

Flood Risk PRA is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Flood Risk PRA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Flood Risk PRA.  If not, see <https://www.gnu.org/licenses/>.

"""

# imports
import os
import shutil
import datetime as dt
import sys
import numpy as np
#import scipy.stats as sstats
import subprocess
import shapely

# parameters
# 3,583 is the maximum realization + flood index count
# 1,000 is the maximum number of realizations.
# 501.3 is the maximum input discharge ...
START_REAL = 51
END_REAL = 100
MAX_REAL = 1000
NROWS = 200
NCOLS = 70
DEPTH_CUTOFF = 0.01
#   directories and file names
IN_PRECIP_DIR = r'C:\Users\nicholas.martin\Documents\Flood_Risk_Model\Fr' \
                r'io_Synthetic_Weather\Processed_Outputs'
IN_PRE_XLSX = "All_Events_Summary-Processed.xlsx"
MOD_FILES_DIR = "Model_Files"
SCRIPT_FILES_DIR = "Py_Scripts"
INPUTS = "input.txt"
DEPTH = "Depth.txt"
TOPO = "Topo.txt"
MANN = "Mann.txt"
CALC_DEPTH = "H.txt"
LOG_FILE = "FR-PRA_Log_R%04dto%04d.txt"
#OUT_MOD = ""
RESULTS_DIR = "Results"
V_FILE = "V.txt"
U_FILE = "U.txt"
#   obstruction location and information
#   commented out for this branch.
#OBS_LOC = ( [167,35], [167,36] )
#OBS_AVAIL_HEIGHT = ( (105.033-91.700), (105.033-91.700) )
#OBS_ORIG_DEPTH = ( 5.01, 5.01 )
# set-up the distribution and sampler for obstruction amount.
#OBS_DEF_SEED = int( 62379 )
#OBS_SAMPLER = None
#OBS_GEV = sstats.genextreme( -0.1, loc=0, scale=0.5 )
# tracking lists
WATER_DEPTH_LIST = list()
FLOOD_DEPTH_LIST = list()
U_VEL_LIST = list()
V_VEL_LIST = list()
# inflow boundary information
INFLOW_TOPO = { 28 : 111.667,
                29 : 110.000,
                30 : 107.500,
                31 : 105.000,
                32 : 103.333,
                33 : 101.667,
                34 : 100.000,
                35 : 100.000,
                36 : 101.667,
                37 : 103.333,
                38 : 105.000,
                39 : 107.500,
                40 : 110.000,
                41 : 111.667, }
INFLOW_BOUND = { 0 : [200.0, 100.0, 107.50, [31, 32, 33, 34, 35, 36, 37, 38,], ],
                 1 : [325.0, 200.0, 110.0, [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, ], ],
                 2 : [425.0, 325.0, 111.667, [29, 30, 31, 32, 33, 34, 35, 36, 37,
                                              38, 39, 40, ], ],
                 3 : [541.7, 425.0, 113.333, [28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, ], ], }
# building information
BUILDING_META = { 0 : [ 1, [ 109.633, ( (112, 27), (113, 27), (114,27), (115,26), ),
                             ( (112,25), (112,26), (113,25), (113,26), (114,26),), ], ],
                  1 : [ 2, [ 109.383, ( (117, 27), (118, 27), (119,27), (120,26), ),
                             ( (117,25), (117,26), (118,25), (118,26), (119,26),), ], ],
                  2 : [ 3, [ 109.133, ( (122, 27), (123, 27), (124,27), (125,26), ),
                             ( (122,25), (122,26), (123,25), (123,26), (124,26),), ], ],
                  3 : [ 4, [ 108.883, ( (127, 27), (128, 27), (129,27), (130,26), ),
                             ( (127,25), (127,26), (128,25), (128,26), (129,26),), ], ],
                  4 : [ 5, [ 108.633, ( (132, 27), (133, 27), (134,27), (135,26), ),
                             ( (132,25), (132,26), (133,25), (133,26), (134,26),), ], ],
                  5 : [ 6, [ 108.383, ( (137, 27), (138, 27), (139,27), (140,26), ),
                             ( (137,25), (137,26), (138,25), (138,26), (139,26),), ], ],
                  6 : [ 7, [ 108.133, ( (142, 27), (143, 27), (144,27), (145,26), ),
                             ( (142,25), (142,26), (143,25), (143,26), (144,26),), ], ],
                  7 : [ 8, [ 107.883, ( (147, 27), (148, 27), (149,27), (150,26), ),
                             ( (147,25), (147,26), (148,25), (148,26), (149,26),), ], ],
                  8 : [ 9, [ 107.633, ( (152, 27), (153, 27), (154,27), (155,26), ),
                             ( (152,25), (152,26), (153,25), (153,26), (154,26),), ], ],
                  9 : [ 10, [ 107.383, ( (157, 27), (158, 27), (159,27), (160,26), ),
                              ( (157,25), (157,26), (158,25), (158,26), (159,26),), ], ],
                 10 : [ 11, [ 107.133, ( (162, 27), (163, 27), (164,27), (165,26), ),
                              ( (162,25), (162,26), (163,25), (163,26), (164,26),), ], ],
                 11 : [ 12, [ 104.30, ( (112, 31), (113, 31), (114,31), (115,30), ),
                              ( (112,29), (112,30), (113,29), (113,30), (114,30),), ], ],
                 12 : [ 13, [ 104.050, ( (117, 31), (118, 31), (119,31), (120,30), ),
                              ( (117,29), (117,30), (118,29), (118,30), (119,30),), ], ],
                 13 : [ 14, [ 103.800, ( (122, 31), (123, 31), (124,31), (125,30), ),
                              ( (122,29), (122,30), (123,29), (123,30), (124,30),), ], ],
                 14 : [ 15, [ 103.550, ( (127, 31), (128, 31), (129,31), (130,30), ),
                              ( (127,29), (127,30), (128,29), (128,30), (129,30),), ], ],
                 15 : [ 16, [ 103.300, ( (132, 31), (133, 31), (134,31), (135,30), ),
                              ( (132,29), (132,30), (133,29), (133,30), (134,30),), ], ],
                 16 : [ 17, [ 103.050, ( (137, 31), (138, 31), (139,31), (140,30), ),
                              ( (137,29), (137,30), (138,29), (138,30), (139,30),), ], ],
                 17 : [ 18, [ 102.800, ( (142, 31), (143, 31), (144,31), (145,30), ),
                              ( (142,29), (142,30), (143,29), (143,30), (144,30),), ], ],
                 18 : [ 19, [ 102.550, ( (147, 31), (148, 31), (149,31), (150,30), ),
                              ( (147,29), (147,30), (148,29), (148,30), (149,30),), ], ],
                 19 : [ 20, [ 102.300, ( (152, 31), (153, 31), (154,31), (155,30), ),
                              ( (152,29), (152,30), (153,29), (153,30), (154,30),), ], ],
                 20 : [ 21, [ 102.050, ( (157, 31), (158, 31), (159,31), (160,30), ),
                              ( (157,29), (157,30), (158,29), (158,30), (159,30),), ], ],
                 21 : [ 22, [ 101.800, ( (162, 31), (163, 31), (164,31), (165,30), ),
                              ( (162,29), (162,30), (163,29), (163,30), (164,30),), ], ],
                 22 : [ 23, [ 104.30, ( (112, 40), (113, 40), (114,40), (115,41), ),
                              ( (112,41), (112,42), (113,41), (113,42), (114,41),), ], ],
                 23 : [ 24, [ 104.050, ( (117, 40), (118, 40), (119,40), (120,41), ),
                              ( (117,41), (117,42), (118,41), (118,42), (119,41),), ], ],
                 24 : [ 25, [ 103.800, ( (122, 40), (123, 40), (124,40), (125,41), ),
                              ( (122,41), (122,42), (123,41), (123,42), (124,41),), ], ],
                 25 : [ 26, [ 103.550, ( (127, 40), (128, 40), (129,40), (130,41), ),
                              ( (127,41), (127,42), (128,41), (128,42), (129,41),), ], ],
                 26 : [ 27, [ 103.300, ( (132, 40), (133, 40), (134,40), (135,41), ),
                              ( (132,41), (132,42), (133,41), (133,42), (134,41),), ], ],
                 27 : [ 28, [ 103.050, ( (137, 40), (138, 40), (139,40), (140,41), ),
                              ( (137,41), (137,42), (138,41), (138,42), (139,41),), ], ],
                 28 : [ 29, [ 102.800, ( (142, 40), (143, 40), (144,40), (145,41), ),
                              ( (142,41), (142,42), (143,41), (143,42), (144,41),), ], ],
                 29 : [ 30, [ 102.550, ( (147, 40), (148, 40), (149,40), (150,41), ),
                              ( (147,41), (147,42), (148,41), (148,42), (149,41),), ], ],
                 30 : [ 31, [ 102.300, ( (152, 40), (153, 40), (154,40), (155,41), ),
                              ( (152,41), (152,42), (153,41), (153,42), (154,41),), ], ],
                 31 : [ 32, [ 102.050, ( (157, 40), (158, 40), (159,40), (160,41), ),
                              ( (157,41), (157,42), (158,41), (158,42), (159,41),), ], ],
                 32 : [ 33, [ 101.800, ( (162, 40), (163, 40), (164,40), (165,41), ),
                              ( (162,41), (162,42), (163,41), (163,42), (164,41),), ], ],
                 33 : [ 34, [ 109.633, ( (112, 44), (113, 44), (114,44), (115,45), ),
                              ( (112,45), (112,46), (113,45), (113,46), (114,45),), ], ],
                 34 : [ 35, [ 109.383, ( (117, 44), (118, 44), (119,44), (120,45), ),
                              ( (117,45), (117,46), (118,45), (118,46), (119,45),), ], ],
                 35 : [ 36, [ 109.133, ( (122, 44), (123, 44), (124,44), (125,45), ),
                              ( (122,45), (122,46), (123,45), (123,46), (124,45),), ], ],
                 36 : [ 37, [ 108.883, ( (127, 44), (128, 44), (129,44), (130,45), ),
                              ( (127,45), (127,46), (128,45), (128,46), (129,45),), ], ],
                 37 : [ 38, [ 108.633, ( (132, 44), (133, 44), (134,44), (135,45), ),
                             ( (132,45), (132,46), (133,45), (133,46), (134,45),), ], ],
                 38 : [ 39, [ 108.383, ( (137, 44), (138, 44), (139,44), (140,45), ),
                              ( (137,45), (137,46), (138,45), (138,46), (139,45),), ], ],
                 39 : [ 40, [ 108.133, ( (142, 44), (143, 44), (144,44), (145,45), ),
                              ( (142,45), (142,46), (143,45), (143,46), (144,45),), ], ],
                 40 : [ 41, [ 107.883, ( (147, 44), (148, 44), (149,44), (150,45), ),
                              ( (147,45), (147,46), (148,45), (148,46), (149,45),), ], ],
                 41 : [ 42, [ 107.633, ( (152, 44), (153, 44), (154,44), (155,45), ),
                              ( (152,45), (152,46), (153,45), (153,46), (154,45),), ], ],
                 42 : [ 43, [ 107.383, ( (157, 44), (158, 44), (159,44), (160,45), ),
                              ( (157,45), (157,46), (158,45), (158,46), (159,45),), ], ],
                 43 : [ 44, [ 107.133, ( (162, 44), (163, 44), (164,44), (165,45), ),
                             ( (162,45), (162,46), (163,45), (163,46), (164,45),), ], ], }
# want to use topo height for cells in valuelist [1][1] and compare that height to valuelist[1][0] to determine if
# have enough water depth for inundation.
NUM_BUILDS = len( BUILDING_META )
BUILDING_POLYS = [ shapely.geometry.Polygon( ( ( 120.0, 555.0 ), ( 130.0, 555.0 ), (130.0, 570.0), (125.0, 570.0 ),
                                             ( 125.0, 565.0 ), ( 120.0, 565.0 ), ( 120.0, 555.0 ), ) ), #0
                   shapely.geometry.Polygon( ( ( 120.0, 580.0 ), ( 130.0, 580.0 ), (130.0, 595.0), (125.0, 595.0 ),
                                             ( 125.0, 590.0 ), ( 120.0, 590.0 ), ( 120.0, 580.0 ), ) ), #1
                   shapely.geometry.Polygon( ( ( 120.0, 605.0 ), ( 130.0, 605.0 ), (130.0, 620.0), (125.0, 620.0 ),
                                             ( 125.0, 615.0 ), ( 120.0, 615.0 ), ( 120.0, 605.0 ), ) ), #2
                   shapely.geometry.Polygon( ( ( 120.0, 630.0 ), ( 130.0, 630.0 ), (130.0, 645.0), (125.0, 645.0 ),
                                             ( 125.0, 640.0 ), ( 120.0, 640.0 ), ( 120.0, 630.0 ), ) ), #3
                   shapely.geometry.Polygon( ( ( 120.0, 655.0 ), ( 130.0, 655.0 ), (130.0, 670.0), (125.0, 670.0 ),
                                             ( 125.0, 665.0 ), ( 120.0, 665.0 ), ( 120.0, 655.0 ), ) ), #4
                   shapely.geometry.Polygon( ( ( 120.0, 680.0 ), ( 130.0, 680.0 ), (130.0, 695.0), (125.0, 695.0 ),
                                             ( 125.0, 690.0 ), ( 120.0, 690.0 ), ( 120.0, 680.0 ), ) ), #5
                   shapely.geometry.Polygon( ( ( 120.0, 705.0 ), ( 130.0, 705.0 ), (130.0, 720.0), (125.0, 720.0 ),
                                             ( 125.0, 715.0 ), ( 120.0, 715.0 ), ( 120.0, 705.0 ), ) ), #6
                   shapely.geometry.Polygon( ( ( 120.0, 730.0 ), ( 130.0, 730.0 ), (130.0, 745.0), (125.0, 745.0 ),
                                             ( 125.0, 740.0 ), ( 120.0, 740.0 ), ( 120.0, 730.0 ), ) ), #7
                   shapely.geometry.Polygon( ( ( 120.0, 755.0 ), ( 130.0, 755.0 ), (130.0, 770.0), (125.0, 770.0 ),
                                             ( 125.0, 765.0 ), ( 120.0, 765.0 ), ( 120.0, 755.0 ), ) ), #8
                   shapely.geometry.Polygon( ( ( 120.0, 780.0 ), ( 130.0, 780.0 ), (130.0, 795.0), (125.0, 795.0 ),
                                             ( 125.0, 790.0 ), ( 120.0, 790.0 ), ( 120.0, 780.0 ), ) ), #9
                   shapely.geometry.Polygon( ( ( 120.0, 805.0 ), ( 130.0, 805.0 ), (130.0, 820.0), (125.0, 820.0 ),
                                             ( 125.0, 815.0 ), ( 120.0, 815.0 ), ( 120.0, 805.0 ), ) ), #10
                   shapely.geometry.Polygon( ( ( 140.0, 555.0 ), ( 150.0, 555.0 ), (150.0, 570.0), (145.0, 570.0 ),
                                             ( 145.0, 565.0 ), ( 140.0, 565.0 ), ( 140.0, 555.0 ), ) ), #11
                   shapely.geometry.Polygon( ( ( 140.0, 580.0 ), ( 150.0, 580.0 ), (150.0, 595.0), (145.0, 595.0 ),
                                             ( 145.0, 590.0 ), ( 140.0, 590.0 ), ( 140.0, 580.0 ), ) ), #12
                   shapely.geometry.Polygon( ( ( 140.0, 605.0 ), ( 150.0, 605.0 ), (150.0, 620.0), (145.0, 620.0 ),
                                             ( 145.0, 615.0 ), ( 140.0, 615.0 ), ( 140.0, 605.0 ), ) ), #13
                   shapely.geometry.Polygon( ( ( 140.0, 630.0 ), ( 150.0, 630.0 ), (150.0, 645.0), (145.0, 645.0 ),
                                             ( 145.0, 640.0 ), ( 140.0, 640.0 ), ( 140.0, 630.0 ), ) ), #14
                   shapely.geometry.Polygon( ( ( 140.0, 655.0 ), ( 150.0, 655.0 ), (150.0, 670.0), (145.0, 670.0 ),
                                             ( 145.0, 665.0 ), ( 140.0, 665.0 ), ( 140.0, 655.0 ), ) ), #15
                   shapely.geometry.Polygon( ( ( 140.0, 680.0 ), ( 150.0, 680.0 ), (150.0, 695.0), (145.0, 695.0 ),
                                             ( 145.0, 690.0 ), ( 140.0, 690.0 ), ( 140.0, 680.0 ), ) ), #16
                   shapely.geometry.Polygon( ( ( 140.0, 705.0 ), ( 150.0, 705.0 ), (150.0, 720.0), (145.0, 720.0 ),
                                             ( 145.0, 715.0 ), ( 140.0, 715.0 ), ( 140.0, 705.0 ), ) ), #17
                   shapely.geometry.Polygon( ( ( 140.0, 730.0 ), ( 150.0, 730.0 ), (150.0, 745.0), (145.0, 745.0 ),
                                             ( 145.0, 740.0 ), ( 140.0, 740.0 ), ( 140.0, 730.0 ), ) ), #18
                   shapely.geometry.Polygon( ( ( 140.0, 755.0 ), ( 150.0, 755.0 ), (150.0, 770.0), (145.0, 770.0 ),
                                             ( 145.0, 765.0 ), ( 140.0, 765.0 ), ( 140.0, 755.0 ), ) ), #19
                   shapely.geometry.Polygon( ( ( 140.0, 780.0 ), ( 150.0, 780.0 ), (150.0, 795.0), (145.0, 795.0 ),
                                             ( 145.0, 790.0 ), ( 140.0, 790.0 ), ( 140.0, 780.0 ), ) ), #20
                   shapely.geometry.Polygon( ( ( 140.0, 805.0 ), ( 150.0, 805.0 ), (150.0, 820.0), (145.0, 820.0 ),
                                             ( 145.0, 815.0 ), ( 140.0, 815.0 ), ( 140.0, 805.0 ), ) ), #21
                   shapely.geometry.Polygon( ( ( 200.0, 555.0 ), ( 210.0, 555.0 ), (210.0, 565.0), (205.0, 565.0 ),
                                             ( 205.0, 570.0 ), ( 200.0, 570.0 ), ( 200.0, 555.0 ), ) ), #22
                   shapely.geometry.Polygon( ( ( 200.0, 580.0 ), ( 210.0, 580.0 ), (210.0, 590.0), (205.0, 590.0 ),
                                             ( 205.0, 595.0 ), ( 200.0, 595.0 ), ( 200.0, 580.0 ), ) ), #23
                   shapely.geometry.Polygon( ( ( 200.0, 605.0 ), ( 210.0, 605.0 ), (210.0, 615.0), (205.0, 615.0 ),
                                             ( 205.0, 620.0 ), ( 200.0, 620.0 ), ( 200.0, 605.0 ), ) ), #24
                   shapely.geometry.Polygon( ( ( 200.0, 630.0 ), ( 210.0, 630.0 ), (210.0, 640.0), (205.0, 640.0 ),
                                             ( 205.0, 645.0 ), ( 200.0, 645.0 ), ( 200.0, 630.0 ), ) ), #25
                   shapely.geometry.Polygon( ( ( 200.0, 655.0 ), ( 210.0, 655.0 ), (210.0, 665.0), (205.0, 665.0 ),
                                             ( 205.0, 670.0 ), ( 200.0, 670.0 ), ( 200.0, 655.0 ), ) ), #26
                   shapely.geometry.Polygon( ( ( 200.0, 680.0 ), ( 210.0, 680.0 ), (210.0, 690.0), (205.0, 690.0 ),
                                             ( 205.0, 695.0 ), ( 200.0, 695.0 ), ( 200.0, 680.0 ), ) ), #27
                   shapely.geometry.Polygon( ( ( 200.0, 705.0 ), ( 210.0, 705.0 ), (210.0, 715.0), (205.0, 715.0 ),
                                             ( 205.0, 720.0 ), ( 200.0, 720.0 ), ( 200.0, 705.0 ), ) ), #28
                   shapely.geometry.Polygon( ( ( 200.0, 730.0 ), ( 210.0, 730.0 ), (210.0, 740.0), (205.0, 740.0 ),
                                             ( 205.0, 745.0 ), ( 200.0, 745.0 ), ( 200.0, 730.0 ), ) ), #29
                   shapely.geometry.Polygon( ( ( 200.0, 755.0 ), ( 210.0, 755.0 ), (210.0, 765.0), (205.0, 765.0 ),
                                             ( 205.0, 770.0 ), ( 200.0, 770.0 ), ( 200.0, 755.0 ), ) ), #30
                   shapely.geometry.Polygon( ( ( 200.0, 780.0 ), ( 210.0, 780.0 ), (210.0, 790.0), (205.0, 790.0 ),
                                             ( 205.0, 795.0 ), ( 200.0, 795.0 ), ( 200.0, 780.0 ), ) ), #31
                   shapely.geometry.Polygon( ( ( 200.0, 805.0 ), ( 210.0, 805.0 ), (210.0, 815.0), (205.0, 815.0 ),
                                             ( 205.0, 820.0 ), ( 200.0, 820.0 ), ( 200.0, 805.0 ), ) ), #32
                   shapely.geometry.Polygon( ( ( 220.0, 555.0 ), ( 230.0, 555.0 ), (230.0, 565.0), (225.0, 565.0 ),
                                             ( 225.0, 570.0 ), ( 220.0, 570.0 ), ( 220.0, 555.0 ), ) ), #33
                   shapely.geometry.Polygon( ( ( 220.0, 580.0 ), ( 230.0, 580.0 ), (230.0, 590.0), (225.0, 590.0 ),
                                             ( 225.0, 595.0 ), ( 220.0, 595.0 ), ( 220.0, 580.0 ), ) ), #34
                   shapely.geometry.Polygon( ( ( 220.0, 605.0 ), ( 230.0, 605.0 ), (230.0, 615.0), (225.0, 615.0 ),
                                             ( 225.0, 620.0 ), ( 220.0, 620.0 ), ( 220.0, 605.0 ), ) ), #35
                   shapely.geometry.Polygon( ( ( 220.0, 630.0 ), ( 230.0, 630.0 ), (230.0, 640.0), (225.0, 640.0 ),
                                             ( 225.0, 645.0 ), ( 220.0, 645.0 ), ( 220.0, 630.0 ), ) ), #36
                   shapely.geometry.Polygon( ( ( 220.0, 655.0 ), ( 230.0, 655.0 ), (230.0, 665.0), (225.0, 665.0 ),
                                             ( 225.0, 670.0 ), ( 220.0, 670.0 ), ( 220.0, 655.0 ), ) ), #37
                   shapely.geometry.Polygon( ( ( 220.0, 680.0 ), ( 230.0, 680.0 ), (230.0, 690.0), (225.0, 690.0 ),
                                             ( 225.0, 695.0 ), ( 220.0, 695.0 ), ( 220.0, 680.0 ), ) ), #38
                   shapely.geometry.Polygon( ( ( 220.0, 705.0 ), ( 230.0, 705.0 ), (230.0, 715.0), (225.0, 715.0 ),
                                             ( 225.0, 720.0 ), ( 220.0, 720.0 ), ( 220.0, 705.0 ), ) ), #39
                   shapely.geometry.Polygon( ( ( 220.0, 730.0 ), ( 230.0, 730.0 ), (230.0, 740.0), (225.0, 740.0 ),
                                             ( 225.0, 745.0 ), ( 220.0, 745.0 ), ( 220.0, 730.0 ), ) ), #40
                   shapely.geometry.Polygon( ( ( 220.0, 755.0 ), ( 230.0, 755.0 ), (230.0, 765.0), (225.0, 765.0 ),
                                             ( 225.0, 770.0 ), ( 220.0, 770.0 ), ( 220.0, 755.0 ), ) ), #41
                   shapely.geometry.Polygon( ( ( 220.0, 780.0 ), ( 230.0, 780.0 ), (230.0, 790.0), (225.0, 790.0 ),
                                             ( 225.0, 795.0 ), ( 220.0, 795.0 ), ( 220.0, 780.0 ), ) ), #42
                   shapely.geometry.Polygon( ( ( 220.0, 805.0 ), ( 230.0, 805.0 ), (230.0, 815.0), (225.0, 815.0 ),
                                             ( 225.0, 820.0 ), ( 220.0, 820.0 ), ( 220.0, 805.0 ), ) ), ] #43


# functions
def readRealizations( LogFile ):
    """Read in the realizations to a DataFrame

    Parameters
    ----------
    LogFile : str
        FQDN log file name.

    Returns
    -------
    RealDF : pd.DataFrame
        Table that was worksheet with events.

    """
    # imports
    import pandas as pd
    # globals
    global IN_PRECIP_DIR, IN_PRE_XLSX
    # parameters
    # locals
    Infiler = os.path.normpath( os.path.join( IN_PRECIP_DIR, IN_PRE_XLSX ) )
    RealDF = pd.read_excel( Infiler, sheet_name="Events", header=0,
                            index_col=0, )
    return RealDF


def adjustInflowBnds( inputFile, curDischarge, LogFile ):
    """Adjust the inflow boundary condition specification to agree with discharge.

    Parameters
    ----------
    inputFile : str
        FQDN name for current inflow.
    curDischarge : float
        Discharge in cms.
    LogFile : str
        FQDN log file name.

    Returns
    -------
    retStatus: 0 == success.

    """
    # imports
    from copy import deepcopy
    # globals
    global INFLOW_TOPO, INFLOW_BOUND, NCOLS
    # parameters
    goodReturn = 0
    badReturn = 1
    KW_IN_VEL = "VELDYVEL"
    KW_IN_DEP = "TDEPDYDEP"
    CELL_LENGTH = 5.0
    # locals
    # start
    # find the current inflow
    NumBndSpec = len( INFLOW_BOUND )
    bndIndex = -1
    for iI in range( NumBndSpec ):
        maxDischarge = INFLOW_BOUND[iI][0]
        minDischarge = INFLOW_BOUND[iI][1]
        if (curDischarge > minDischarge) and (curDischarge <= maxDischarge):
            # then have found it
            bndIndex = iI
            break
        # end if
    # end for
    # check
    if bndIndex == -1:
        # Then this is an error
        errMsg = "Did not find boundary specification for discharge %6.2f!!!\n" % curDischarge
        with open( LogFile, 'a' ) as LF:
            LF.write( "%s" % errMsg )
        # end if
        return badReturn
    # end if
    # if made it here then need to read and process input file to set new
    #   inflow boundaries.
    with open( inputFile, 'r' ) as Inf:
        AllLines = Inf.readlines()
    # end with
    OutLines = deepcopy( AllLines )
    # now process
    lCnt = 0
    for tLine in AllLines:
        stripLine = tLine.strip()
        if len(stripLine) < 3:
            lCnt += 1
            continue
        if stripLine[0] == "#":
            lCnt += 1
            continue
        # end checks
        if "=" in stripLine:
            # see if need to process.
            initSplit = stripLine.split("=")
            if initSplit[0].strip() == KW_IN_DEP:
                newDepths = [ 0.0 for x in range(NCOLS) ]
                setLocs = INFLOW_BOUND[bndIndex][3]
                numLocs = len( setLocs )
                lowerDischarge = INFLOW_BOUND[bndIndex][1]
                refElevation = INFLOW_BOUND[bndIndex][2]
                topHeight = ( ( curDischarge - lowerDischarge) /
                             ( float( numLocs ) * CELL_LENGTH * 1.0 ) )
                adjElev = INFLOW_TOPO[setLocs[0]]
                for cC in range(numLocs):
                    if ( cC == 0 ) or ( cC == numLocs-1 ):
                        newDepths[setLocs[cC]] = topHeight
                    else:
                        newDepths[setLocs[cC]]  = ( ( adjElev -
                                                      INFLOW_TOPO[setLocs[cC]] )
                                                    + topHeight )
                    # end if
                # end for
                # make the new line
                newLiner = "%s = " % KW_IN_DEP
                for cC in range(NCOLS):
                    newLiner += "%5.2f " % newDepths[cC]
                # end for
                newLiner += "\n"
                # now assign
                OutLines[lCnt] = newLiner
                lCnt += 1
            elif initSplit[0] == KW_IN_VEL:
                newVels = [ 0.0 for x in range(NCOLS) ]
                setLocs = INFLOW_BOUND[bndIndex][3]
                for cC in setLocs:
                    newVels[cC] = 1.0
                # end for
                # make the new line
                newLiner = "%s = " % KW_IN_VEL
                for cC in range(NCOLS):
                    newLiner += "%5.2f " % newVels[cC]
                # end for
                newLiner += "\n"
                # now assign
                OutLines[lCnt] = newLiner
                lCnt += 1
            else:
                lCnt += 1
                continue
            # end if
        else:
            lCnt += 1
            continue
        # end if
    # end for
    # now output the modified
    with open( inputFile, 'w' ) as Inf:
        Inf.writelines(OutLines)
    # end with
    # return
    return goodReturn


def adjustDepthandTopo( depFile, topoFile, curObs, LogFile ):
    """Adjust the water depth for the obstruction


    Parameters
    ----------
    depFile : str
        String name for active depth file.
    topoFile : str
        String name for current topography file
    curObs : float
        Current obstruction depth as sampled from OBS_GEV
    LogFile : str
        String name for log file.

    Returns
    -------
    retStatus : int
        0 == success, anything else is failure

    """
    # imports
    # globals
    global OBS_LOC, OBS_AVAIL_HEIGHT, OBS_ORIG_DEPTH, NROWS, NCOLS
    # parameters
    goodReturn = 0
    badReturn = -1
    # locals
    # start
    # calculate the available depth
    availDepth = OBS_AVAIL_HEIGHT[0]
    if curObs > availDepth:
        useObs = availDepth
    else:
        useObs = curObs
    # end if
    newTopo = np.loadtxt( topoFile )
    newH = np.loadtxt( depFile )
    # check
    if ( newTopo.shape[0] != NROWS ) or ( newH.shape[0] != NROWS ):
        with open( LogFile, 'a' ) as LF:
            LF.write( "Either topo or depth file was loaded incorrectly!!!")
        # end with
        return badReturn
    # end if
    if ( newTopo.shape[1] != NCOLS ) or ( newH.shape[1] != NCOLS ):
        with open( LogFile, 'a' ) as LF:
            LF.write( "Either topo or depth file was loaded incorrectly!!!")
        # end with
        return badReturn
    # end if
    # if here then can process
    existDepth1 = newH[OBS_LOC[0][0]-1, OBS_LOC[0][1]-1]
    newDepth1 = existDepth1 - useObs
    if newDepth1 <= 0.0:
        newDepth1 = 0.0
    # end if
    newH[OBS_LOC[0][0]-1, OBS_LOC[0][1]-1] = newDepth1
    existDepth2 = newH[OBS_LOC[1][0]-1, OBS_LOC[1][1]-1]
    newDepth2 = existDepth2 - useObs
    if newDepth2 <= 0.0:
        newDepth2 = 0.0
    # end if
    newH[OBS_LOC[1][0]-1, OBS_LOC[1][1]-1] = newDepth2
    existTopo1 = newTopo[OBS_LOC[0][0]-1, OBS_LOC[0][1]-1]
    newTopo1 = existTopo1 + useObs
    newTopo[OBS_LOC[0][0]-1, OBS_LOC[0][1]-1] = newTopo1
    existTopo2 = newTopo[OBS_LOC[1][0]-1, OBS_LOC[1][1]-1]
    newTopo2 = existTopo2 + useObs
    newTopo[OBS_LOC[1][0]-1, OBS_LOC[1][1]-1] = newTopo2
    # now write out depth
    with open(depFile, 'w+') as OF:
        for iI in range(NROWS):
            for jJ in range(NCOLS):
                OF.write('%6.2f   ' % newH[iI,jJ])
            # end for
            OF.write("\n")
        # end for
    # end with
    with open(topoFile, 'w+') as OF:
        for iI in range(NROWS):
            for jJ in range(NCOLS):
                OF.write('%6.2f   ' % newTopo[iI,jJ])
            # end for
            OF.write("\n")
        # end for
    # end with
    # return
    return goodReturn


def processFlooding( CWD, realNum, floodNum, curObs, curDis, LogFile ):
    """Determine flooding for this realization

    Parameters
    ----------
    CWD : str
        Current working directory, where files are
    realNum : int
        Realization number or index.
    floodNum : int
        Flood number or index within this realization
    curObs : float
        Obstruction depth.
    curDis : float
        Current input discharge
    LogFile : str
        Log file name.

    Returns
    -------
    inunDF : pd.DataFrame

    """
    # imports
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    # globals
    global TOPO, CALC_DEPTH, DEPTH_CUTOFF, NROWS, NCOLS, RESULTS_DIR
    global NUM_BUILDS, BUILDING_META, BUILDING_POLYS, V_FILE, U_FILE
    global FLOOD_DEPTH_LIST, U_VEL_LIST, V_VEL_LIST, WATER_DEPTH_LIST
    # parameters
    # locals
    # start
    cTopoFile = os.path.normpath( os.path.join( CWD, TOPO ) )
    topo = np.loadtxt( cTopoFile, dtype=np.float32, )
    cDepFile = os.path.normpath( os.path.join( CWD, CALC_DEPTH ) )
    H1Array = np.loadtxt( cDepFile, dtype=np.float32 )
    H1Array = np.where( H1Array <= DEPTH_CUTOFF, 0.0, H1Array )
    H = np.reshape( H1Array, (NROWS, NCOLS), order='C' ).copy()
    # calculate inundation
    InBuildDict = dict()
    for cB in range(NUM_BUILDS):
        cBId = BUILDING_META[cB][0]
        cFoundElev = BUILDING_META[cB][1][0]
        checkLocTuple = BUILDING_META[cB][1][1][2]
        checkLocTopo = topo[checkLocTuple[0]-1, checkLocTuple[1]-1]
        cFoundHeight = cFoundElev - checkLocTopo
        cWaterDepth = H[checkLocTuple[0]-1, checkLocTuple[1]-1]
        cInunDepth = cWaterDepth - cFoundHeight
        if cInunDepth < 0.0:
            cInunDepth = 0.0
        # end if
        # add to our output dictionary
        InBuildDict[cBId] = [ checkLocTuple[0], checkLocTuple[1],
                             checkLocTopo, cFoundElev, cFoundHeight,
                             cWaterDepth, cInunDepth ]
    # end for
    # make a DataFrame
    DataDict = { "Row" : np.array( [ InBuildDict[x][0] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.int32 ),
                 "Column" : np.array( [ InBuildDict[x][1] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.int32 ),
                 "Topo_m" : np.array( [ InBuildDict[x][2] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.float32 ),
                 "FloorEl_m" : np.array( [ InBuildDict[x][3] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.float32 ),
                 "FloorHeight_m" : np.array( [ InBuildDict[x][4] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.float32 ),
                 "WaterDepth_m" : np.array( [ InBuildDict[x][5] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.float32 ),
                 "FloodDepth_m" : np.array( [ InBuildDict[x][6] for x in range(1, NUM_BUILDS+1, 1) ], dtype=np.float32 ), }
    InunDF = pd.DataFrame( index=[x for x in range(1, NUM_BUILDS+1, 1)], data=DataDict )
    # make a plot
    # make mesh grid
    InFiler = os.path.normpath( os.path.join( CWD, "XINDEX.txt" ) )
    XINDEX = np.loadtxt( InFiler, dtype=np.float32, )
    InFiler = os.path.normpath( os.path.join( CWD, "YINDEX.txt" ) )
    YINDEX = np.loadtxt( InFiler, dtype=np.float32, )
    X_Pts = np.array( [ XINDEX[i] + ( 0.5*(XINDEX[i+1] - XINDEX[i]) )
                        for i in range(NCOLS) ], dtype=np.float32)
    Y_Pts = np.array( [ YINDEX[j] + ( 0.5*(YINDEX[j+1] - YINDEX[j]) )
                        for j in range(NROWS) ], dtype=np.float32)
    YV, XV = np.meshgrid( Y_Pts, X_Pts, indexing='ij' )
    # need to get U and V
    InFiler = os.path.normpath( os.path.join( CWD, U_FILE ) )
    U1Array = np.loadtxt( InFiler, dtype=np.float32, )
    InFiler = os.path.normpath( os.path.join( CWD, "Hux.txt" ) )
    Hux1Array = np.loadtxt( InFiler, dtype=np.float32, )
    Hux1Array = np.where( Hux1Array <= DEPTH_CUTOFF, 0.0, Hux1Array )
    U1Array = np.where( Hux1Array <= DEPTH_CUTOFF, 0.0, U1Array )
    InFiler = os.path.normpath( os.path.join( CWD, V_FILE ) )
    V1Array = np.loadtxt( InFiler, dtype=np.float32, )
    InFiler = os.path.normpath( os.path.join( CWD, "Hvy.txt" ) )
    Hvy1Array = np.loadtxt( InFiler, dtype=np.float32, )
    Hvy1Array = np.where( Hvy1Array <= DEPTH_CUTOFF, 0.0, Hvy1Array )
    V1Array = np.where( Hvy1Array <= DEPTH_CUTOFF, 0.0, V1Array )
    # put on regular grid
    npU = np.zeros( (NROWS*NCOLS), dtype=np.float32 )
    npV = np.zeros( (NROWS*NCOLS), dtype=np.float32 )
    for iI in range( NROWS ):
        for jJ in range(NCOLS):
            cUm1 = ( ( iI * (NCOLS+1) ) + jJ )
            cUp1 = ( ( iI * (NCOLS+1) ) + ( jJ+1 ) )
            cVm1 = ( ( iI * NCOLS ) + jJ )
            cVp1 = ( ( ( iI + 1 ) * NCOLS ) + jJ )
            curNode = ( ( iI * NCOLS ) + jJ )
            npU[curNode] = 0.5 * ( U1Array[cUm1] + U1Array[cUp1] )
            npV[curNode] = 0.5 * ( V1Array[cVm1] + V1Array[cVp1] )
        # end for
    # end for
    npU = np.where( np.abs( npU ) < 0.0001, 0.0, npU )
    npV = np.where( np.abs( npV ) < 0.0001, 0.0, npV )
    plotU = np.reshape( npU, (NROWS, NCOLS), order='C' ).copy()
    plotV = np.reshape( npV, (NROWS, NCOLS), order='C' ).copy()
    # contour intervals for water depth
    CLevels = np.array( [x*0.5 for x in range(41)], dtype=np.float32 )
    CLevels[0] = 0.01
    # plot coordinates
    xplotR = [ 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0 ]
    yplotR = [ 400.0, 500.0, 600.0, 700.0, 800.0, 900.0 ]
    # output file name
    OutFiler = "R%04d_Fl%02d_Focus_Area_WLVel.png" % ( realNum, floodNum )
    OutFilePNG = os.path.normpath( os.path.join( CWD, RESULTS_DIR, OutFiler ) )
    # plot
    Fig1 = plt.figure()
    Fig1.set_size_inches(5.0, 8.0)
    ax11 = Fig1.add_subplot(1,1,1)
    wd = ax11.contourf( XV, YV, H, levels=CLevels, cmap='Blues', vmin=0.01, vmax=20.0, zorder=20.0)
    qv = ax11.quiver( XV[::3, 28:42:1], YV[::3, 28:42:1], plotU[::3, 28:42:1], plotV[::3, 28:42:1],
                      headwidth=5, color='xkcd:tangerine', scale=30, zorder=50.0 )
    ax11.set_ylabel( "Northing (m)", fontsize=10)
    ax11.set_xlabel( "Easting (m)", fontsize=10)
    ax11.set_xticks( xplotR )
    ax11.set_yticks( yplotR )
    cb = Fig1.colorbar( wd, ax=ax11, orientation='vertical',)
    cb.ax.tick_params(labelsize=8)
    cb.set_label( 'Water Depth (m)', fontsize=9 )
    pCnt = 0
    for tPoly in BUILDING_POLYS:
        zLev = pCnt + 25.0
        zAnno = pCnt + 50.0
        pX, pY = tPoly.exterior.xy
        cCCoords = tPoly.representative_point().coords[:][0]
        labelstr = "%d" % (pCnt+1)
        ax11.fill( pX, pY, edgecolor='xkcd:medium grey', facecolor='xkcd:medium grey', linewidth=1, zorder=zLev, alpha=1.0 )
        ax11.text( cCCoords[0]-5.0, cCCoords[1]-5.0, labelstr, color='xkcd:cyan', fontsize=7, fontweight='normal', zorder=zAnno )
        pCnt += 1
    # end for
    ax11.grid( visible=True, which='major', axis='y' )
    ax11.set_xlim( (xplotR[0], xplotR[len(xplotR)-1]) )
    ax11.set_ylim( (yplotR[0], yplotR[len(yplotR)-1]) )
    ax11.tick_params(axis='both', which='major', labelsize=9)
    ax11.xaxis.set_major_formatter( mpl.ticker.StrMethodFormatter( "{x:,.0f}" ) )
    ax11.yaxis.set_major_formatter( mpl.ticker.StrMethodFormatter( "{x:,.0f}" ) )
    Fig1.savefig( OutFilePNG, dpi=600 )
    Fig1.clf()
    plt.close(fig=Fig1)
    # assign to our collation/summary tracking lists
    WATER_DEPTH_LIST.append( float( InunDF["WaterDepth_m"].max() ) )
    FLOOD_DEPTH_LIST.append( float( InunDF["FloodDepth_m"].max() ) )
    U_VEL_LIST.append( float( npU.max() ) )
    V_VEL_LIST.append( float( npV.max() ) )
    # return
    return InunDF


def outputSummary( CWD, ClRealList, FlIndList, DTList, PrecipList, DisList,
                   ObsDepList, FloodDFList, LogFile ):
    """Output the inundation and input configuration summary for these realizations.

    Parameters
    ----------
    CWD : str
        current working directory.
    ClRealList : list
        List of integers for climate realizations.
    FlIndList : list
        List of integers for flood indexes
    DTList : list
        List of pd.Timestamps for projected date of the flood.
    PrecipList : list
        List of floats for precipitation depth.
    DisList : list
        List of float discharges for scaled input discharge.
    ObsDepList : list
        List of floats for obstruction depth
    FloodDFList : list of pd.DataFrame
        List of DataFrames with the inundation summary
    LogFile : str
        FQDN for log file.

    Returns
    -------
    None.

    """
    # imports
    import pandas as pd
    # globals
    global FLOOD_DEPTH_LIST, U_VEL_LIST, V_VEL_LIST, WATER_DEPTH_LIST
    global START_REAL, END_REAL
    # globals
    # parameters
    # locals
    # start
    OutFiler = "R%04dto%04d_Flooding_Summary_All.xlsx" % (START_REAL, END_REAL )
    OutFP = os.path.normpath( os.path.join( CWD, RESULTS_DIR, OutFiler ) )
    DataDict = { "Realization" : np.array( ClRealList, dtype=np.int32 ),
                 "Flood Num." : np.array( FlIndList, dtype=np.int32 ),
                 "Date" : DTList,
                 "Precip_mm" : np.array( PrecipList, dtype=np.float32 ),
                 "Discharge_cms" : np.array( DisList, dtype=np.float32 ),
                 "Obstruction_Depth_m" : np.array( ObsDepList, dtype=np.float32 ),
                 "Max_Water_Depth_m" : np.array( WATER_DEPTH_LIST, dtype=np.float32 ),
                 "Max_Flood_Depth_m" : np.array( FLOOD_DEPTH_LIST, dtype=np.float32 ),
                 "Max_U_mps" : np.array( U_VEL_LIST, dtype=np.float32 ),
                 "Max_V_mps" : np.array( V_VEL_LIST, dtype=np.float32 ),}
    SummaryDF = pd.DataFrame( data=DataDict )
    # output to Excel
    writer = pd.ExcelWriter( OutFP )
    workbook  = writer.book
    format1 = workbook.add_format({'num_format': '#,##0.000'})
    format2 = workbook.add_format({'num_format': '#,##0.0000'})
    format3 = workbook.add_format({'num_format': '#,##0'})
    cLabel = "Summary"
    SummaryDF.to_excel( writer, sheet_name=cLabel, )
    # adjust columns
    writer.sheets[cLabel].set_column( 0, 0, 18 )
    for column in SummaryDF:
        column_width = max(SummaryDF[column].astype(str).map(len).max()+6, len(column)+6)
        col_idx = SummaryDF.columns.get_loc(column)
        if column in ["Realization", "Flood Num."]:
            writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width, format3)
        elif column in ["Date"]:
            writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width, )
        else:
            writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width, format1)
        # end if
    # end for
    # now output all of our inundation calculations
    NumFlReal = len( FloodDFList )
    for iI in range(NumFlReal):
        cLabel = "Inun_R%04d_Fl%02d" % ( ClRealList[iI], FlIndList[iI] )
        curDF = FloodDFList[iI]
        curDF.to_excel( writer, sheet_name=cLabel, )
        # adjust columns
        writer.sheets[cLabel].set_column( 0, 0, 10 )
        for column in curDF:
            column_width = max(curDF[column].astype(str).map(len).max()+6, len(column)+6)
            col_idx = curDF.columns.get_loc(column)
            if column in ["Row", "Column"]:
                writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width, format3)
            elif column in ["Topo_m", "FloorEl_m", "FloorHeight_m", "WaterDepth_m", "FloodDepth_m", ]:
                writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width, format1)
            else:
                writer.sheets[cLabel].set_column(col_idx+1, col_idx+1, column_width,)
            # end if
        # end for
    # end for real flood pairs
    writer.close()
    # return
    return


#standalone execution block
# assumes that this module is executed within the same current directory
# as the input file
if __name__ == "__main__":
    CWD = os.getcwd()
    LogFile = os.path.normpath( os.path.join( CWD, LOG_FILE % ( START_REAL, END_REAL ) ) )
    StartDT = dt.datetime.now()
    with open( LogFile, 'w+' ) as LF:
        LF.write( "Start of Flood Risk, Probabilistic Risk Assessment (PRA) - no blockages \n")
        LF.write( "Start time: %s \n\n" % StartDT.strftime("%Y-%m-%d %H:%M"))
    # end with
    MFilesDir = os.path.normpath( os.path.join( CWD, MOD_FILES_DIR ) )
    RealDF = readRealizations( LogFile )
    # initialize tracking structures.
    ClRealList = list()
    FlIndList = list()
    DTList = list()
    PrecipList = list()
    DisList = list()
    ObsDepList = list()
    FloodDFList = list()
    # Now do by number of realizations
    for rR in range(START_REAL, END_REAL+1):
        # get the climate realization and use to set the seed and random sampler
        #curSeed = OBS_DEF_SEED + rR
        #OBS_SAMPLER = np.random.RandomState( seed=curSeed )
        # get floods for only this realization
        curRealDF = RealDF[RealDF["RealNum"] == rR].copy()
        # check to make sure that there are floods
        if len(curRealDF) <= 0:
            # log message
            OutStr = "Climate realization %d has 0 floods.\n" % rR
            with open( LogFile, 'a' ) as LF:
                LF.write("%s" % OutStr )
            # end with
            # then continue
            continue
        # end if
        # if made it here then have floods.
        flCnt = 1
        for indx, row in curRealDF.iterrows():
            curInDischarge = float( row["Discharge_cms"] )
            ClRealList.append( rR )
            FlIndList.append( flCnt )
            DTList.append( row["DateTime"] )
            PrecipList.append( float( row["Precip_mm"] ) )
            DisList.append( curInDischarge )
            # copy base files
            srcFile = os.path.normpath( os.path.join( MFilesDir, INPUTS ) )
            newInFile = os.path.normpath( os.path.join( CWD, INPUTS ) )
            try:
                oF = shutil.copyfile( srcFile, newInFile )
            except:
                OutStr = "Error copying file %s to %s !!!\n" % ( srcFile, newInFile )
                with open( LogFile, 'a' ) as LF:
                    LF.write("%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end try
            srcFile = os.path.normpath( os.path.join( MFilesDir, DEPTH ) )
            newDepFile = os.path.normpath( os.path.join( CWD, DEPTH ) )
            try:
                oF = shutil.copyfile( srcFile, newDepFile )
            except:
                OutStr = "Error copying file %s !!!\n" % DEPTH
                with open( LogFile, 'a' ) as LF:
                    LF.write("%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end try
            srcFile = os.path.normpath( os.path.join( MFilesDir, TOPO ) )
            newTopoFile = os.path.normpath( os.path.join( CWD, TOPO ) )
            try:
                oF = shutil.copyfile( srcFile, newTopoFile )
            except:
                OutStr = "Error copying file %s !!!\n" % TOPO
                with open( LogFile, 'a' ) as LF:
                    LF.write("%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end try
            srcFile = os.path.normpath( os.path.join( MFilesDir, MANN ) )
            dstFile = os.path.normpath( os.path.join( CWD, MANN ) )
            try:
                oF = shutil.copyfile( srcFile, dstFile )
            except:
                OutStr = "Error copying file %s !!!\n" % MANN
                with open( LogFile, 'a' ) as LF:
                    LF.write("%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end try
            # update the input file for the new discharge.
            retStatus = adjustInflowBnds( newInFile, curInDischarge, LogFile )
            if retStatus != 0:
                # then there was an error
                with open( LogFile, 'a' ) as LF:
                    OutStr = "Error in realization %d writing inflow boundary!!!\n" % rR
                    LF.write( "%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end if
            # get the current obstruction depth and set the seed for this realization
            curObstruction = 0.0
            #curObstruction = float( OBS_GEV.rvs( size=1,
            #                                     random_state=OBS_SAMPLER )[0] )
            #if curObstruction < 0.0:
            #    curObstruction = 0.0
            # end if
            ObsDepList.append( curObstruction )
            # write entry to the log file
            with open( LogFile, 'a' ) as LF:
                LF.write( "Climate realization %d, flood index %d, obstruction " \
                          "depth %5.2f, discharge %6.2f \n" %
                          (rR, flCnt, curObstruction, curInDischarge) )
            # end with
            # modify the depth file to reflect the obstruction
            #retStatus = adjustDepthandTopo( newDepFile, newTopoFile,
            #                                curObstruction, LogFile )
            #if retStatus != 0:
            #    # then there was an error
            #    with open( LogFile, 'a' ) as LF:
            #        OutStr = "Error in realization %d writing updated topo" \
            #                 " and depth!!!\n" % rR
            #        LF.write( "%s" % OutStr )
            #    # end with
            #    sys.exit([-1, OutStr])
            ## end if
            # now run
            runResult = subprocess.run( ["MOD_FreeSurf2D.exe"], shell=True,
                                        capture_output=True, text=True, )
            if runResult.returncode != 0:
                # then there was an error
                with open( LogFile, 'a' ) as LF:
                    LF.write( "%s\n\n" % print(runResult.stdout) )
                    LF.write( "%s\n\n" % print(runResult.stderr) )
                # end with
                sys.exit([-1, "Error in MOD_FreeSurf2D execution"])
            # end if
            # process results
            curFloodDF = processFlooding( CWD, rR, flCnt, curObstruction,
                                          curInDischarge, LogFile )
            if len( curFloodDF ) <= 0:
                # then there was an error
                with open( LogFile, 'a' ) as LF:
                    OutStr = "Error in climate realization %d, flood index %d " \
                             "collating outputs!!!\n" % (rR, flCnt)
                    LF.write( "%s" % OutStr )
                # end with
                sys.exit([-1, OutStr])
            # end if
            # add to the tracking list
            FloodDFList.append( curFloodDF )
            # increment the counter
            flCnt += 1
        # end of flood index for
    # end of climate realization for
    # output summary info
    outputSummary( CWD, ClRealList, FlIndList, DTList, PrecipList, DisList,
                   ObsDepList, FloodDFList, LogFile )
    # log file wrap up
    EndDT = dt.datetime.now()
    ETimeDelta = EndDT - StartDT
    ETimeHrs = ( ETimeDelta.total_seconds() / (60.0*60.0) )
    with open( LogFile, 'a' ) as LF:
        OutStr = "Successful completion at %s, elapsed time %6.2f hours \n" % (
                     EndDT.strftime("%Y-%m-%d %H:%M"), ETimeHrs )
        LF.write( "%s" % OutStr )
    # end with
    # done

#EOF