import time
import math
import numpy as np
import os
import shutil
from multiprocessing import Process
from tempfile import mkdtemp
from GribFile import GribFile
from thermo import rh_to_mr
import ozone as o3

_101_pressure_levels = \
                   [0.0050,    0.0161,    0.0384,    0.0769,    0.1370, \
         0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972, \
         1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204, \
         5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492, \
        14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829, \
        29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882, \
        51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396, \
        83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775, \
       125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784, \
       180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338, \
       247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369, \
       328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738, \
       424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200, \
       535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398, \
       661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897, \
       802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236, \
       958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000]

_101_pressure_values = [0] * len(_101_pressure_levels)

#for i in range(101):
    #_101_pressure_levels[i] = [_101_pressure_levels[i], -1]



grb = None
nx = 9
ny = 35
nz = 3
init = True
med_cached = False
coef = np.zeros((nx + 1, ny, nz))
cc   = np.zeros((nz, nx + 1, ny))
tempPres = None
rhPres = None
outPres = None
coordGrid = None
NUM_PROCS = 64

#  Zone 1

#  0.0050 mb:
coef[:, 0,0] = [768.137268,   -0.294463,   -1.938120,    1.445014,    1.010484, \
                 -1.501080,   -0.333852,   -0.460916,   -0.289109,   -0.185883]
#  0.0161 mb:
coef[:, 1,0] = [729.598755,   -0.048072,   -2.192048,    1.556883,    1.041545, \
                 -1.527906,   -0.179924,   -0.397472,   -0.527085,   -0.022882]
#  0.0384 mb:
coef[:, 2,0] = [580.007996,    0.256773,   -1.953103,    1.262224,    0.855626, \
                 -1.226672,    0.021504,   -0.272799,   -0.672050,    0.167226]
#  0.0769 mb:
coef[:, 3,0] = [378.418610,    0.481383,   -1.547368,    0.949718,    0.557970, \
                 -0.797499,    0.236666,   -0.090441,   -0.724859,    0.324352]
#  0.1370 mb:
coef[:, 4,0] = [227.726013,    0.399807,   -0.767517,    0.339597,    0.334890, \
                 -0.304244,    0.242581,    0.200534,   -0.710240,    0.359780]
#  0.2244 mb:
coef[:, 5,0] = [131.540237,    0.116343,    0.181285,   -0.404714,    0.194162, \
                  0.168128,    0.112673,    0.505991,   -0.649249,    0.319506]
#  0.3454 mb:
coef[:, 6,0] = [36.948456,   -0.219132,    1.068293,   -1.006580,    0.041164, \
                 0.596757,    0.016176,    0.724429,   -0.524707,    0.280223]
#  0.5064 mb:
coef[:, 7,0] = [-47.064129,   -0.502069,    1.832718,   -1.521932,   -0.102700, \
                  0.973025,   -0.065917,    0.916755,   -0.417642,    0.247873]
#  0.7140 mb:
coef[:, 8,0] = [-125.058113,   -0.371252,    1.929798,   -1.497271,   -0.445293, \
                   1.206479,   -0.047493,    1.051199,   -0.411198,    0.284401]
#  0.9753 mb:
coef[:, 9,0] = [-195.826279,   -0.252554,    2.017885,   -1.474895,   -0.756146, \
                   1.418304,   -0.030777,    1.173188,   -0.405352,    0.317545]
#  1.2972 mb:
coef[:,10,0] = [-211.840836,   -0.087453,    1.892584,   -1.498378,   -0.888045, \
                   1.557249,   -0.013888,    1.153322,   -0.274422,    0.228628]
#  1.6872 mb:
coef[:,11,0] = [-222.487411,    0.077679,    1.743851,   -1.517234,   -0.992018, \
                   1.676924,    0.003067,    1.119941,   -0.144600,    0.140581]
#  2.1526 mb:
coef[:,12,0] = [-220.345123,    0.284957,    1.478203,   -1.409153,   -1.098716, \
                   1.732799,    0.005315,    1.062679,   -0.051710,    0.083439]
#  2.7009 mb:
coef[:,13,0] = [-191.945847,    0.575167,    0.992410,   -1.054349,   -1.230327, \
                   1.674808,   -0.024700,    0.961055,   -0.019512,    0.071966]
#  3.3398 mb:
coef[:,14,0] = [-163.447739,    0.789576,    0.662758,   -0.848763,   -1.165136, \
                   1.486376,   -0.005209,    0.856173,    0.036504,    0.000516]
#  4.0770 mb:
coef[:,15,0] = [-135.140701,    0.938526,    0.466544,   -0.770585,   -0.932123, \
                   1.188257,    0.055934,    0.750415,    0.112469,   -0.122608]
#  4.9204 mb:
coef[:,16,0] = [-110.534798,    1.079643,    0.269122,   -0.685459,   -0.724919, \
                   0.927309,    0.105440,    0.666502,    0.180002,   -0.240317]
#  5.8776 mb:
coef[:,17,0] = [-86.496407,    1.197924,    0.180591,   -0.680802,   -0.560434, \
                  0.725634,    0.114811,    0.524036,    0.254776,   -0.292023]
#  6.9567 mb:
coef[:,18,0] = [-63.631668,    1.308646,    0.105890,   -0.683525,   -0.407437, \
                  0.538658,    0.120171,    0.383015,    0.326700,   -0.335437]
#  8.1655 mb:
coef[:,19,0] = [-53.883522,    1.446655,   -0.084576,   -0.584522,   -0.321630, \
                  0.418490,    0.138251,    0.323507,    0.266197,   -0.297677]
#  9.5119 mb:
coef[:,20,0] = [-45.056870,    1.579389,   -0.270613,   -0.486306,   -0.242174, \
                  0.306220,    0.155973,    0.269677,    0.203612,   -0.258671]
# 11.0038 mb:
coef[:,21,0] = [-37.774498,    1.605633,   -0.336989,   -0.423904,   -0.179642, \
                  0.233162,    0.156346,    0.226620,    0.162179,   -0.225690]
# 12.6492 mb:
coef[:,22,0] = [-31.381693,    1.580467,   -0.344828,   -0.379910,   -0.126498, \
                  0.180357,    0.148424,    0.189608,    0.131713,   -0.196274]
# 14.4559 mb:
coef[:,23,0] = [-25.256893,    1.556356,   -0.352337,   -0.337762,   -0.075582, \
                  0.129767,    0.140834,    0.154147,    0.102526,   -0.168092]
# 16.4318 mb:
coef[:,24,0] = [-19.332832,    1.532327,   -0.358523,   -0.296751,   -0.023514, \
                  0.076359,    0.135129,    0.118184,    0.075598,   -0.141048]
# 18.5847 mb:
coef[:,25,0] = [-13.621389,    1.508887,   -0.364070,   -0.257120,    0.027774, \
                  0.023140,    0.130262,    0.082869,    0.050141,   -0.115058]
# 20.9224 mb:
coef[:,26,0] = [-10.552977,    1.475920,   -0.355611,   -0.209939,    0.047681, \
                 -0.006736,    0.119757,    0.072153,    0.024709,   -0.098524]
# 23.4526 mb:
coef[:,27,0] = [-11.404002,    1.427825,   -0.325822,   -0.150293,    0.020678, \
                 -0.002060,    0.100506,    0.098321,   -0.001260,   -0.095886]
# 26.1829 mb:
coef[:,28,0] = [-12.431881,    1.380402,   -0.295441,   -0.093749,   -0.004390, \
                  0.002642,    0.081415,    0.123695,   -0.026677,   -0.092560]
# 29.1210 mb:
coef[:,29,0] = [-13.700622,    1.333231,   -0.263907,   -0.040468,   -0.027290, \
                  0.007437,    0.062283,    0.148370,   -0.051712,   -0.088305]
# 32.2744 mb:
coef[:,30,0] = [-12.049180,    1.288188,   -0.233483,   -0.020884,   -0.028961, \
                  0.007528,    0.048806,    0.133084,   -0.050300,   -0.074679]
# 35.6505 mb:
coef[:,31,0] = [-9.317835,    1.244821,   -0.204065,   -0.014504,   -0.022521, \
                 0.005827,    0.037740,    0.102881,   -0.038848,   -0.057748]
# 39.2566 mb:
coef[:,32,0] = [-6.672400,    1.202820,   -0.175573,   -0.008325,   -0.016283, \
                 0.004179,    0.027023,    0.073628,   -0.027756,   -0.041350]
# 43.1001 mb:
coef[:,33,0] = [-4.107985,    1.162105,   -0.147954,   -0.002335,   -0.010237, \
                 0.002583,    0.016633,    0.045271,   -0.017004,   -0.025454]
# 47.1882 mb:
coef[:,34,0] = [-1.620239,    1.122606,   -0.121160,    0.003476,   -0.004369, \
                 0.001033,    0.006554,    0.017761,   -0.006572,   -0.010032]
# Zone 2

#  0.0050 mb:
coef[:, 0,1] = [834.861145,   -0.395518,   -0.168292,   -0.522989,    0.135368, \
                 -0.648098,   -0.120294,   -0.432757,   -0.616573,   -0.049551]
#  0.0161 mb:
coef[:, 1,1] = [706.697998,   -0.564473,    0.061239,   -0.490616,    0.279786, \
                 -0.515674,   -0.074882,   -0.145787,   -0.887622,    0.150642]
#  0.0384 mb:
coef[:, 2,1] = [456.640137,   -0.515665,    0.238243,   -0.325252,    0.306180, \
                 -0.235547,   -0.028916,    0.203108,   -0.970324,    0.303589]
#  0.0769 mb:
coef[:, 3,1] = [196.650543,   -0.580368,    0.469164,   -0.167309,    0.389823, \
                  0.028119,    0.028069,    0.512286,   -0.967577,    0.455387]
#  0.1370 mb:
coef[:, 4,1] = [55.473385,   -0.577106,    0.510232,   -0.037092,    0.481069, \
                 0.095268,    0.059780,    0.602277,   -0.792924,    0.480604]
#  0.2244 mb:
coef[:, 5,1] = [3.620781,   -0.523341,    0.420839,    0.076334,    0.556607, \
                0.033327,    0.063618,    0.596428,   -0.589134,    0.438588]
#  0.3454 mb:
coef[:, 6,1] = [-45.342590,   -0.570448,    0.387072,    0.196123,    0.643398, \
                 -0.025933,    0.053534,    0.719213,   -0.558915,    0.459968]
#  0.5064 mb:
coef[:, 7,1] = [-88.503006,   -0.608382,    0.353316,    0.301524,    0.715639, \
                 -0.073664,    0.042024,    0.832244,   -0.535806,    0.479966]
#  0.7140 mb:
coef[:, 8,1] = [-120.005302,   -0.541255,    0.223263,    0.373417,    0.655807, \
                   0.010448,   -0.035588,    1.041446,   -0.612080,    0.524942]
#  0.9753 mb:
coef[:, 9,1] = [-148.589035,   -0.480346,    0.105258,    0.438649,    0.601519, \
                   0.086768,   -0.106009,    1.231265,   -0.681287,    0.565751]
#  1.2972 mb:
coef[:,10,1] = [-161.973740,   -0.371849,    0.008552,    0.486533,    0.557450, \
                   0.112235,   -0.175632,    1.253408,   -0.627670,    0.561722]
#  1.6872 mb:
coef[:,11,1] = [-173.125641,   -0.266075,   -0.083040,    0.532244,    0.515987, \
                   0.132649,   -0.241772,    1.260636,   -0.568181,    0.555640]
#  2.1526 mb:
coef[:,12,1] = [-173.081696,   -0.094634,   -0.220959,    0.567489,    0.477990, \
                   0.120870,   -0.283719,    1.244042,   -0.533282,    0.546141]
#  2.7009 mb:
coef[:,13,1] = [-150.789734,    0.220479,   -0.455112,    0.578471,    0.446871, \
                   0.041562,   -0.277308,    1.177748,   -0.543711,    0.525773]
#  3.3398 mb:
coef[:,14,1] = [-130.673721,    0.435536,   -0.553035,    0.536432,    0.402999, \
                  -0.011230,   -0.249925,    1.061322,   -0.498618,    0.491557]
#  4.0770 mb:
coef[:,15,1] = [-112.541351,    0.564179,   -0.533735,    0.448692,    0.348486, \
                  -0.040873,   -0.204408,    0.902095,   -0.405841,    0.445350]
#  4.9204 mb:
coef[:,16,1] = [-96.168854,    0.685395,   -0.516674,    0.364636,    0.299450, \
                 -0.066547,   -0.160329,    0.753376,   -0.318513,    0.400629]
#  5.8776 mb:
coef[:,17,1] = [-121.107452,    0.789144,   -0.466781,    0.308981,    0.265531, \
                  -0.041420,   -0.096691,    0.627374,   -0.172447,    0.330783]
#  6.9567 mb:
coef[:,18,1] = [-148.557831,    0.886490,   -0.416297,    0.258454,    0.234542, \
                  -0.012946,   -0.034285,    0.509281,   -0.027978,    0.261965]
#  8.1655 mb:
coef[:,19,1] = [-146.946228,    0.924881,   -0.351371,    0.173360,    0.243558, \
                  -0.023024,   -0.043987,    0.511167,   -0.028737,    0.240058]
#  9.5119 mb:
coef[:,20,1] = [-144.347504,    0.959377,   -0.288869,    0.090871,    0.253623, \
                  -0.034051,   -0.055879,    0.517345,   -0.034760,    0.220859]
# 11.0038 mb:
coef[:,21,1] = [-136.725281,    0.967216,   -0.250049,    0.038883,    0.268749, \
                  -0.045579,   -0.034804,    0.465370,   -0.010830,    0.193404]
# 12.6492 mb:
coef[:,22,1] = [-126.862816,    0.962160,   -0.223352,    0.002548,    0.285976, \
                  -0.057106,    0.001579,    0.386704,    0.026907,    0.162579]
# 14.4559 mb:
coef[:,23,1] = [-117.413834,    0.957316,   -0.197774,   -0.032263,    0.302480, \
                  -0.068150,    0.036436,    0.311337,    0.063062,    0.133046]
# 16.4318 mb:
coef[:,24,1] = [-108.387199,    0.951026,   -0.170067,   -0.067915,    0.319109, \
                  -0.079884,    0.071843,    0.243068,    0.098536,    0.099678]
# 18.5847 mb:
coef[:,25,1] = [-99.728241,    0.944342,   -0.142208,   -0.103052,    0.335396, \
                 -0.091603,    0.106632,    0.179040,    0.132930,    0.065653]
# 20.9224 mb:
coef[:,26,1] = [-90.187675,    0.953354,   -0.131290,   -0.116374,    0.326715, \
                 -0.094256,    0.128964,    0.132450,    0.141525,    0.044614]
# 23.4526 mb:
coef[:,27,1] = [-79.100502,    0.986260,   -0.145695,   -0.097073,    0.280154, \
                 -0.083285,    0.132998,    0.111127,    0.111379,    0.042697]
# 26.1829 mb:
coef[:,28,1] = [-68.699837,    1.017853,   -0.157891,   -0.083327,    0.241075, \
                 -0.074152,    0.137479,    0.089334,    0.082811,    0.041165]
# 29.1210 mb:
coef[:,29,1] = [-59.048897,    1.048161,   -0.167399,   -0.076558,    0.211126, \
                 -0.067266,    0.142594,    0.066655,    0.055907,    0.040111]
# 32.2744 mb:
coef[:,30,1] = [-48.262180,    1.062744,   -0.160901,   -0.062914,    0.173613, \
                 -0.056035,    0.123367,    0.051704,    0.041447,    0.034117]
# 35.6505 mb:
coef[:,31,1] = [-37.250641,    1.071060,   -0.148434,   -0.046917,    0.133942, \
                 -0.043364,    0.095243,    0.039982,    0.032004,    0.026359]
# 39.2566 mb:
coef[:,32,1] = [-26.585524,    1.079114,   -0.136360,   -0.031422,    0.095519, \
                 -0.031093,    0.068005,    0.028629,    0.022858,    0.018845]
# 43.1001 mb:
coef[:,33,1] = [-16.247211,    1.086921,   -0.124655,   -0.016402,    0.058273, \
                 -0.019197,    0.041600,    0.017624,    0.013991,    0.011561]
# 47.1882 mb:
coef[:,34,1] = [-6.217585,    1.094496,   -0.113301,   -0.001831,    0.022139, \
                -0.007657,    0.015985,    0.006947,    0.005390,    0.004494]

# Zone 3

#  0.0050 mb:
coef[:, 0,2] = [253.099335,   -0.059435,   -0.113692,   -0.076584,    0.018085, \
                  0.109696,   -0.245047,    0.207798,   -0.158821,   -0.005658]
#  0.0161 mb:
coef[:, 1,2] = [393.724976,   -0.159516,   -0.305137,   -0.205541,    0.048536, \
                  0.294412,   -0.657678,    0.557706,   -0.426257,   -0.015186]
#  0.0384 mb:
coef[:, 2,2] = [488.316467,   -0.220589,   -0.421964,   -0.284236,    0.067118, \
                  0.407132,   -0.909481,    0.771234,   -0.589456,   -0.021001]
#  0.0769 mb:
coef[:, 3,2] = [571.788940,   -0.274628,   -0.525335,   -0.353867,    0.083561, \
                  0.506871,   -1.132283,    0.960168,   -0.733859,   -0.026146]
#  0.1370 mb:
coef[:, 4,2] = [564.131775,   -0.280373,   -0.502551,   -0.299793,    0.065189, \
                  0.444018,   -1.027926,    0.888957,   -0.690901,   -0.011248]
#  0.2244 mb:
coef[:, 5,2] = [496.362030,   -0.249809,   -0.394284,   -0.179865,    0.027459, \
                  0.290057,   -0.717002,    0.634090,   -0.503432,    0.004950]
#  0.3454 mb:
coef[:, 6,2] = [432.359314,   -0.211118,   -0.286078,   -0.105599,   -0.000052, \
                  0.184573,   -0.442425,    0.356257,   -0.266806,   -0.011249]
#  0.5064 mb:
coef[:, 7,2] = [375.968506,   -0.173903,   -0.194453,   -0.041546,   -0.020229, \
                  0.090298,   -0.207287,    0.122826,   -0.062737,   -0.026678]
#  0.7140 mb:
coef[:, 8,2] = [335.498596,   -0.064613,   -0.226756,   -0.032027,    0.072557, \
                 -0.012652,   -0.217699,    0.255347,   -0.032882,   -0.068283]
#  0.9753 mb:
coef[:, 9,2] = [298.778290,    0.034553,   -0.256067,   -0.023391,    0.156749, \
                 -0.106065,   -0.227148,    0.375591,   -0.005793,   -0.106034]
#  1.2972 mb:
coef[:,10,2] = [320.508484,    0.055737,   -0.215899,   -0.126142,    0.126488, \
                 -0.119170,   -0.163850,    0.312299,    0.057370,   -0.183508]
#  1.6872 mb:
coef[:,11,2] = [338.109192,    0.077756,   -0.167569,   -0.228182,    0.081726, \
                 -0.109376,   -0.114042,    0.249524,    0.124284,   -0.253638]
#  2.1526 mb:
coef[:,12,2] = [330.158142,    0.125628,   -0.127452,   -0.284897,    0.018188, \
                  -0.073173,   -0.058458,    0.164222,    0.202521,   -0.283837]
#  2.7009 mb:
coef[:,13,2] = [288.525238,    0.208034,   -0.113437,   -0.262148,   -0.070341, \
                 -0.019104,    0.050349,   -0.000575,    0.297252,   -0.249445]
#  3.3398 mb:
coef[:,14,2] = [247.402054,    0.288580,   -0.105926,   -0.234689,   -0.143902, \
                  0.044019,    0.059252,   -0.048267,    0.325860,   -0.189348]
#  4.0770 mb:
coef[:,15,2] = [207.174408,    0.365898,   -0.103597,   -0.202744,   -0.204937, \
                  0.116567,   -0.020110,    0.001814,    0.300148,   -0.106910]
#  4.9204 mb:
coef[:,16,2] = [172.771362,    0.425351,   -0.097642,   -0.168375,   -0.266588, \
                  0.200371,   -0.115398,    0.021846,    0.299176,   -0.026388]
#  5.8776 mb:
coef[:,17,2] = [170.819962,    0.377608,    0.003744,   -0.188331,   -0.277131, \
                  0.162590,   -0.022318,   -0.064583,    0.305203,   -0.023369]
#  6.9567 mb:
coef[:,18,2] = [171.848145,    0.322555,    0.108887,   -0.212190,   -0.282633, \
                  0.115752,    0.083180,   -0.156449,    0.311570,   -0.027389]
#  8.1655 mb:
coef[:,19,2] = [178.900009,    0.317759,    0.072287,   -0.147714,   -0.288088, \
                  0.079424,    0.109033,   -0.207073,    0.293082,   -0.003934]
#  9.5119 mb:
coef[:,20,2] = [185.850937,    0.315015,    0.032178,   -0.082944,   -0.293293, \
                  0.045129,    0.130805,   -0.253892,    0.274527,    0.019458]
# 11.0038 mb:
coef[:,21,2] = [181.033875,    0.342789,   -0.008198,   -0.050784,   -0.276437, \
                  0.028043,    0.139765,   -0.267162,    0.249229,    0.043197]
# 12.6492 mb:
coef[:,22,2] = [170.696304,    0.384559,   -0.047857,   -0.034871,   -0.249395, \
                  0.019533,    0.142417,   -0.264131,    0.221237,    0.066608]
# 14.4559 mb:
coef[:,23,2] = [160.792282,    0.424578,   -0.085853,   -0.019626,   -0.223487, \
                  0.011380,    0.144958,   -0.261227,    0.194418,    0.089037]
# 16.4318 mb:
coef[:,24,2] = [150.499008,    0.461093,   -0.122477,   -0.003129,   -0.196627, \
                  0.000506,    0.147210,   -0.251144,    0.167289,    0.109366]
# 18.5847 mb:
coef[:,25,2] = [140.299698,    0.495449,   -0.157735,    0.013451,   -0.170036, \
                 -0.011132,    0.149301,   -0.238611,    0.140674,    0.128437]
# 20.9224 mb:
coef[:,26,2] = [133.652267,    0.528523,   -0.174538,    0.028105,   -0.154952, \
                 -0.021779,    0.150992,   -0.227335,    0.124328,    0.123706]
# 23.4526 mb:
coef[:,27,2] = [132.214737,    0.560411,   -0.163866,    0.040182,   -0.156891, \
                 -0.031171,    0.152119,   -0.217697,    0.123112,    0.082945]
# 26.1829 mb:
coef[:,28,2] = [130.282120,    0.591975,   -0.157004,    0.055184,   -0.159140, \
                 -0.038202,    0.157241,   -0.213196,    0.124561,    0.042131]
# 29.1210 mb:
coef[:,29,2] = [127.686752,    0.623526,   -0.154959,    0.074146,   -0.161817, \
                 -0.042282,    0.167575,   -0.215254,    0.129460,    0.000725]
# 32.2744 mb:
coef[:,30,2] = [108.781189,    0.698895,   -0.147134,    0.068980,   -0.139383, \
                 -0.037117,    0.145935,   -0.184685,    0.111908,   -0.009194]
# 35.6505 mb:
coef[:,31,2] = [84.030518,    0.789494,   -0.137261,    0.054729,   -0.107821, \
                -0.028531,    0.112540,   -0.142284,    0.086147,   -0.006936]
# 39.2566 mb:
coef[:,32,2] = [60.058582,    0.877243,   -0.127698,    0.040927,   -0.077253, \
                -0.020216,    0.080196,   -0.101216,    0.061195,   -0.004749]
# 43.1001 mb:
coef[:,33,2] = [36.821239,    0.962303,   -0.118428,    0.027547,   -0.047621, \
                -0.012155,    0.048844,   -0.061407,    0.037009,   -0.002629]
# 47.1882 mb:
coef[:,34,2] = [14.277436,    1.044823,   -0.109435,    0.014568,   -0.018874, \
                -0.004335,    0.018428,   -0.022788,    0.013544,   -0.000572]


def _get_lower_bound(vertProf, presVec, pressure):
    '''Given a sorted vertical profile array and a pressure, return the highest
    pressure level in the profile array which is less than the given pressure.
    If no pressure level provides a lower bound, return None.
    '''
    for i in xrange(len(presVec) - 1, -1, -1):
        # find the first pressure level less or equal to than the given 
        # pressure and then return it. this only works when vertProf is sorted
        if presVec[i] <= pressure:
            return (presVec[i], vertProf[i])
    
    return None

def _get_upper_bound(vertProf, presVec, pressure):
    '''Given a sorted vertical profile array and a pressure, return the lowest
    pressure level in the profile array which is greater than the given pressure.
    If no pressure level provides an upper bound, return None.
    '''
    for i in xrange(len(presVec)):
        # find the first pressure level greater than or equal to the given 
        # pressure and then return it. this only works when vertProf is sorted
        if presVec[i] >= pressure:
            return (presVec[i], vertProf[i])
    
    return None
            

def _vert_interp_low_pressure(vertProf, presVec, pressure, lat):
    '''Use regression with parabolic fit coefficients to extrapolate to lower
    pressure levels than are given in the input profile. Modeled after extem101.f
    '''
    global init, med_cached
    x1 = 75.0
    x2 = 45.0
    x3 = 15.0
    lx = [35, 39, 44, 50, 55, 63, 69, 75, 85]
    cfl = np.zeros(nx+1)

    # if we haven't already cached the medium pressure levels, cache them
    if not med_cached:
        # calculate values at indices from lx of the 101-level profile
        for i in xrange(len(lx)):
            _101_pressure_values[lx[i]] = _vert_interp_med_pressure(vertProf, presVec, _101_pressure_levels[lx[i]])

        med_cached = True

    # calculate parabolic fit coefficients
    if init:
        init = False
        for j in xrange(ny):
            for i in xrange(nx + 1):

                y1 = coef[i,j,0]
                y2 = coef[i,j,1]
                y3 = coef[i,j,2]
                
                x12 = x1*x1
                x22 = x2*x2
                x32 = x3*x3

                t1 = x2*x32 - x3*x22
                t2 = -(x1*x32 - x3*x12)
                t3 = x1*x22 - x2*x12

                det = t1 + t2 + t3

                cc[0,i,j] = (y1*t1 + y2*t2 + y3*t3) / det
                cc[1,i,j] = ((y2*x32 - y3*x22) - (y1*x32 - y3*x12) + (y1*x22 - y2*x12)) / det
                cc[2,i,j] = ((x2*y3 - x3*y2) - (x1*y3 - x3*y1) + (x1*y2 - x2*y1)) / det

    # load predictors
    plat = np.zeros(nz)
    alat = abs(lat)
    xx = 1
    for k in xrange(nz):
        plat[k] = xx
        xx = xx * alat

    tx = np.zeros(nx)
    for i in xrange(nx):
        tx[i] = _101_pressure_values[lx[i]]

    for j in xrange(ny):
        for i in xrange(nx + 1):
            cy = 0
            for k in xrange(nz):
                cy = cy + cc[k,i,j] * plat[k]
            cfl[i] = cy
        sum = cfl[0]
        for i in xrange(1, nx+1):
            sum += cfl[i] * tx[i-1]
        _101_pressure_values[j] = sum

    return _vert_interp_med_pressure(_101_pressure_values, _101_pressure_levels, pressure, True)


def _vert_interp_high_pressure(vertProf, presVec, pressure):
    '''Use logarithmic extrapolation using the highest possible pressure estimate
    to get values for pressures higher than those given in the profile.
    '''
    p1 = _get_lower_bound(vertProf, presVec, pressure)
    p0 = _get_lower_bound(vertProf, presVec, p1[0] - 0.001) # get next lowest level

    if pressure == p1[0]:
        return p1[1]

    # extrapolate highest pressure in 101-level pressure profile
    y2 = p1[1] + ((p1[1] - p0[1]) * \
            math.log(float(_101_pressure_levels[-1]) / float(p1[0])) / \
            math.log(float(p1[0]) / float(p0[0])))

    # put extrapolated point in a tuple for interpolation
    p2 = (_101_pressure_levels[-1], y2)

    interpVal = _log_interp(p1, p2, pressure)

    #print 'lower  val: %f, %f' % (p1[0], p1[1])
    #print 'interp val: %f, %f' % (pressure, interpVal)
    #print 'upper  val: %f, %f' % (p2[0], p2[1])

    return interpVal

def _vert_interp_med_pressure(vertProf, presVec, pressure, verbose=False):
    '''Use logarithmic interpolation to calculate values between pressure
    levels in the input profile.
    '''
    p2 = _get_upper_bound(vertProf, presVec, pressure)
    p1 = _get_lower_bound(vertProf, presVec, pressure)

    if pressure == p1[0]:
        return p1[1]
    if pressure == p2[0]:
        return p2[1]

    interpVal = _log_interp(p1, p2, pressure)

    #if verbose:
        #print 'lower  val: %f, %f' % (p1[0], p1[1])
        #print 'interp val: %f, %f' % (pressure, interpVal)
        #print 'upper  val: %f, %f' % (p2[0], p2[1])

    return interpVal


def _log_interp(p1, p2, x):
    '''Perform the actual logarithmic interpolation.'''
    return p1[1] + (((p2[1] - p1[1]) * \
            math.log(float(x) / float(p1[0]))) / \
            math.log(float(p2[0]) / float(p1[0])))


def _get_input_profs(varName, lat, lon, filename, tempProf, rhProf, grid=False):
    global grb

    # need temp prof for both mixing ratio and temperature interpolation
    retDict = {}
    coordGrid = None
    rhPres = None
    tempPres = None
    if tempProf == None:
        # if no profile or file was provided, there has been an error
        if filename == None:
            print '_get_input_profs: no input data provided'
            return None

        # only create a new grib file if it is different from the current one
        if grb is None or grb.filename != filename:
            grb = GribFile(filename)

        if not grid:
            tempPres, tempProf = grb.readDataLatLon('Temperature', lat, lon)
            retDict['tempPres'] = tempPres
            retDict['tempProf'] = tempProf
        else:
            coordGrid, tempPres, tempProf = grb.readDataAllLatLon('Temperature')
            retDict['coordGrid'] = coordGrid
            retDict['tempProf'] = tempProf
            retDict['tempPres'] = tempPres
    else:
        retDict['tempPres'] = tempProf[0]
        retDict['tempProf'] = tempProf[1]

    # if we weren't supplied a rel hum profile and that is what we are trying
    # to measure, get one from the provided file (which should exist)
    if rhProf == None and varName == 'Relative humidity':
        # if no profile or file was provided, there has been an error
        if filename == None:
            print '_get_input_profs: no input data provided'
            return None

        # only create a new grib file if it is different from the current one
        # or none currently exists
        if grb is None or grb.filename != filename:
            grb = GribFile(filename)

        if not grid:
            rhPres, rhProf = grb.readDataLatLon(varName, lat, lon)
            retDict['rhPres'] = rhPres
            retDict['rhProf'] = rhProf
        else:
            coordGrid, rhPres, rhProf = grb.readDataAllLatLon(varName)
            retDict['coordGrid'] = coordGrid
            retDict['rhPres'] = rhPres
            retDict['rhProf'] = rhProf
    elif varName == 'Relative humidity':
        # If rhProf is not none, and we are measuring relative humidity
        retDict['rhPres'] = rhProf[0]
        retDict['rhProf'] = rhProf[1]
    else:
        # We are not measuring relative humidity
        retDict['rhPres'] = None
        retDict['rhProf'] = None


    # We already have a temperature and relative humidity profile, but need a grid.
    # Grab the grid from our input profile
    if tempProf is not None and rhProf is not None and grid:
        if filename == None:
            print '_get_input_profs: no input data provided'
            return None

        # only create a new grib file if it is different from the current one
        if grb is None or grb.filename != filename:
            grb = GribFile(filename)

        # Get the coordinate grid
        coordGrid, tempPres, tempProf = grb.readDataAllLatLon('Temperature')
        retDict['coordGrid'] = coordGrid

    return retDict


def _make_ret_prof(pressure, grid=False, gridSize=None):
    # Copy pressure values into our return profile
    retProf = []
    if pressure == None:
        # if no pressure given, just assume 101-level profile
        retProf = [ [_101_pressure_levels[i], _101_pressure_values[i]] for i in xrange(len(_101_pressure_levels)) ]
        if grid:
            # if our profile is a grid, make the output columns a grid of size (nLat, nLons)
            for i in xrange(len(retProf)):
                retProf[i][1] = np.ndarray((gridSize[0], gridSize[1]))

    elif type(pressure) == int or type(pressure) == float or type(pressure) == np.float64:
        # if the pressure is a single number, convert it to a list for processing
        retProf.append([pressure, -1])
        if grid:
            for i in xrange(len(retProf)):
                retProf[i][1] = np.ndarray((gridSize[0], gridSize[1]))

    elif type(pressure) == list or type(pressure) == np.ndarray:
        # if the pressure was given as a list, copy it into our return 
        # profile for processing
        for p in pressure:
            if type(p) != int and type(p) != float and type(p) != np.float64:
                # make sure the pressure levels are of correct type
                print 'Error: incompatible pressure type \'%s\'' % type(p)
                return None
            retProf.append([p, -1])
        if grid:
            for i in xrange(len(retProf)):
                retProf[i][1] = np.ndarray((gridSize[0], gridSize[1]))
    else:
        # pressure is not a compatible data type
        print 'Error: pressure is an incompatible data type \'%s\'' % type(pressure)
        return None

    return retProf


def vert_interp(varName, lat, lon, pressure=None, filename=None, tempProf=None, rhProf=None):
    '''Given an input profile or file from which to extract one, as well as a latitude,
    longitude, and (list of) pressure level(s), compute either the interpolated temperature 
    or mixing ratio, based on the varName given.

    varName  - the name of the measurement we are looking to interpolate
    lat      - the latitude of the measurement
    lon      - the longitude of the measurement
    pressure - the output pressure vector
    filename - the file we want to read input profiles from
    tempProf - A pair of lists (pressure list, data list), which specifies an 
               input temperature profile
    rhProf   - A pair of lists (pressure list, data list), which specifies an 
               input relative humidity profile

    return: list of [pressure, measurement] pairs
    '''
    global med_cached

    # make us recalculate the medium pressure values the first time we do low-pressure interpolation
    med_cached = False

    # Copy pressure values into our return profile
    retProf = _make_ret_prof(pressure)

    inputDict = _get_input_profs(varName, lat, lon, filename, tempProf, rhProf) 
    tempPres = inputDict['tempPres']
    rhPres = inputDict['rhPres']
    tempProf = inputDict['tempProf']
    rhProf = inputDict['rhProf']

    # temp prof should always have something in it now
    if tempProf is None:
        return None

    if varName == 'Relative humidity':
        # If we are computing relative humidity, first convert to mixing ratio
        # before interpolating
        vertProf = []
        vertProf[:] = rhProf[:]
        for i in xrange(len(vertProf)):
            rh = rhProf[i]
            pres = rhPres[i]
            temp = None

            # make sure we get the temperature for the correct pressure level
            for j in xrange(len(tempPres)):
                if tempPres[j] == pres:
                    temp = tempProf[j]
                    break
            if temp == None:
                # If, for some reason, the temp and rh profiles don't match up, error out
                print 'Error: could not find temperature for pressure level %d' % rhPres[i]
                return None

            # Convert rel hum to mixing ratio
            vertProf[i] = max(rh_to_mr(rh, np.array([pres]), np.array([temp]))[0], 0.0030)

        vertPres = rhPres

    elif varName == 'Temperature':
        # If we are computing Temperature, just grab the temperature vector and
        # interpolate
        vertProf = tempProf
        vertPres = tempPres

    else:
        print 'Error: Unsupported varName: %s' % varName
        exit(1)

    for i in xrange(len(retProf)):
        if retProf[i][0] < vertPres[0]:
            # if the pressure is atmospherically higher than the original 26 levels
            #print 'low pressure:'
            if varName == 'Temperature':
                # use regression to calculate temperature at low pressure
                retProf[i][1] = _vert_interp_low_pressure(tempProf, tempPres, retProf[i][0], lat)
            elif varName == 'Relative humidity':
                # fix the water level at low pressure
                retProf[i][1] = 0.003

        elif retProf[i][0] > vertPres[-1]:
            # if the pressure is atmospherically lower than the original 26 levels
            # extrapolate using the highest possible pressure level
            #print 'high pressure:'
            temp = _vert_interp_high_pressure(tempProf, tempPres, retProf[i][0])
            if varName == 'Temperature':
                retProf[i][1] = temp
            elif varName == 'Relative humidity':
                mr = _vert_interp_high_pressure(vertProf, vertPres, retProf[i][0])
                #print temp, mr
                retProf[i][1] = mr

        else:
            # if the pressure falls within the original 26 levels
            # just interpolate
            #print 'med pressure:'
            temp = _vert_interp_med_pressure(tempProf, tempPres, retProf[i][0], True)
            if varName == 'Temperature':
                retProf[i][1] = temp
            elif varName == 'Relative humidity':
                mr = _vert_interp_med_pressure(vertProf, vertPres, retProf[i][0], True)
                #print temp, mr
                retProf[i][1] = mr

    return retProf

def vert_interp_ozone(filename, lat, lon, month):
    '''Given an input profile extracted using the filename, latitude,
    and longitude, compute the extended ozone profile for
    a location in the given month.
    '''
    global grb

    if grb is None or grb.filename != filename:
        grb = GribFile(filename)

    o3prof = grb.readDataLatLon('O3MR', lat, lon)

    # convert ozone units
    for i in xrange(len(o3prof)):
        o3prof[i] = (o3prof[i][0], o3prof[i][1] * 1000000.0 * (28.97 / 48.0))

    # get estimated profile
    ozone = o3.clozo101(lat, month)
    # adjust profile based on GDAS data
    adj_ozone = o3.adjo3(ozone, o3prof)

    retProf = []
    for i in xrange(len(adj_ozone)):
        retProf.append([_101_pressure_levels[i], float(adj_ozone[i])])

    return retProf

def _process_columns(varName, tempFn, rhFn, outFn, tempRecSz, rhRecSz, \
                     outRecSz, offset, nLats, nLons, nRecs):
    global tempPres
    global rhPres
    global outPres
    global coordGrid

    if varName != 'Temperature' and varName != 'Relative humidity':
        print 'Invalid varName %s' % varName
        exit(1)

    # always get the temperture vectors
    tempFp = np.memmap(tempFn, dtype='float64', mode='c', shape=(nLats * nLons, tempRecSz))

    # only get the relative humidity vectors if we are measuring relative humidity
    rhFp = None
    if varName == 'Relative humidity':
        rhFp = np.memmap(rhFn, dtype='float64', mode='c', shape=(nLats * nLons, rhRecSz))

    outFp = np.memmap(outFn, dtype='float64', shape=(nLats * nLons, outRecSz))

    # process nRecs records
    col = 0
    for i in xrange(offset, offset + nRecs):

            if i >= nLats * nLons:
                # if we don't have any more records to process, 
                # just break out of the loop
                break

            lat, lon = divmod(i, nLons)
            #print lat, lon

            #t = time.time()
            tempData = tempFp[i]
            tempCol = (tempPres, tempData)

            rhCol = None
            if rhFp is not None:
                rhData = rhFp[i]
                rhCol = (rhPres, rhData)

            column = vert_interp(varName, coordGrid[lat, lon, 0], coordGrid[lat, lon, 1], 
                                 outPres, None, tempCol, rhCol)

            #column = [ (0, 0) for i in range(len(retProf)) ]

            for j in xrange(len(column)):
                column[j] = column[j][1]
                #retProf[i][1][lat][lon] = column[i][1]

            #for j in xrange(len(column)):
                #outFp[i, j] = column[j]
            outFp[i, :] = column[:]

            #print outFp[i]

            #print 'prof col %5d: %1.7f' % (col, time.time() - t) 
            col += 1

    #print '-' * 80
    #print outFp[offset]
    #print '-' * 80
    #outFp.flush()
    #del outFp
    return


def vert_interp_grid(varName, pressure=None, filename=None, tempProf=None, rhProf=None):
    global coordGrid
    global tempPres
    global rhPres
    global outPres
    inputDict = _get_input_profs(varName, 0, 0, filename, tempProf, rhProf, True)
    coordGrid = inputDict['coordGrid']
    tempPres = inputDict['tempPres']
    rhPres = inputDict['rhPres']
    tempProf = inputDict['tempProf']
    rhProf = inputDict['rhProf']

    if pressure is None:
        outPres = [ _101_pressure_levels[i] for i in xrange(len(_101_pressure_levels)) ]
    else:
        outPres = pressure

    nLats = len(coordGrid)
    nLons = len(coordGrid[0])

    # reshape our input profiles to make them easier to process
    tempProf = tempProf.reshape((len(tempProf), nLats * nLons)).T
    if rhProf is not None:
        rhProf   = rhProf.reshape((len(rhProf), nLats * nLons)).T

    # temp prof should always have something in it now
    if tempProf is None:
        return None

    gridSize = (nLats, nLons)

    # TODO: Might not need to make a return profile (we can just get the length of outPres)
    # - Our return profile is just outFp
    retProf = _make_ret_prof(outPres, True, gridSize)

    # Can now memmap our input and output profiles
    tmpDir = mkdtemp()

    tempFile = os.path.join(tmpDir, 'temp_file.dat')
    tempFp = np.memmap(tempFile, dtype='float64', mode='w+', shape=tempProf.shape)
    tempFp[:] = tempProf[:]
    del tempFp

    rhFile = None
    if rhProf is not None:
        rhFile = os.path.join(tmpDir, 'rh_file.dat')
        rhFp = np.memmap(rhFile, dtype='float64', mode='w+', shape=rhProf.shape)
        rhFp[:] = rhProf[:]
        del rhFp

    outFile = os.path.join(tmpDir, 'out_file.dat')
    outFp = np.memmap(outFile, dtype='float64', mode='w+', shape=(nLats * nLons, len(retProf)))
    del outFp

    columnsPerProc = ((nLats * nLons) // NUM_PROCS) + 1

    # process columns in parallel
    procList = []
    if rhProf is None:
        rhProfLen = None
    else:
        rhProfLen = len(rhProf[0])
    for i in range(NUM_PROCS): 
        p = Process(target=_process_columns, args=(varName, tempFile, rhFile,
                                                   outFile, len(tempProf[0]),
                                                   rhProfLen, len(retProf), 
                                                   i * columnsPerProc, 
                                                   nLats, nLons, columnsPerProc))
        procList.append(p)
        p.start()

    # wait for running child processes to complete
    for p in procList:
        p.join()
        if p.is_alive():
            print 'Error: process %d should be dead here' % p.pid

    # get the output profile from our mem-mapped file
    dataGrid = np.memmap(outFile, dtype='float64', mode='r', shape=(nLats * nLons, len(retProf)))

    # remove tmpDir when done
    shutil.rmtree(tmpDir)

    dataGrid = dataGrid.T.reshape((len(retProf), nLats, nLons))

    # TODO: might not need to return coordGrid. Just added this for data validation purposes
    # Return (lat x lon) coord grid, (pres) pressure column, and (pres x lat x lon) data grid
    return coordGrid, outPres, dataGrid
