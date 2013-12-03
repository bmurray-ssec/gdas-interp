#!/usr/bin/env python

'''
Graciously supplied by Geoff Cureton

https://svn.ssec.wisc.edu/repos/geoffc/Python/Science/thermo.py

Edited to work on numpy arrays.
'''

from scipy import log10
import numpy as np
import math


def rh_to_mr(rh, p, t):
  '''
  Returns mixing ratio, in g/kg, given relative humidity in %,
  pressure in hPa and temperature in K.
  '''
  return rh * 0.01 * satmix(p, t)

def rh_to_mr_wat( rh, p, t) :
  '''
  Returns mixing ratio over water, in g/kg, given relative humidity in %, 
  pressure in hPa and temperature in K.
  '''
  return rh * 0.01 * satmixwat(p, t)

def rh_to_mr_ice( rh, p, t) :
  '''
  Returns mixing ratio over ice, in g/kg, given relative humidity in %, 
  pressure in hPa and temperature in K.
  '''
  return rh * 0.01 * satmixice(p, t)

def mr_to_rh( mr,  p,  t) :
  '''
  Returns relative humidity in %, given the mixing ratio in g/kg,  
  pressure in hPa and temperature in K.
  '''
  return mr * 100. / satmix(p, t)

def mr_to_rh_wat( mr,  p,  t) :
  '''
  Returns relative humidity in %, given the mixing ratio over water in g/kg,  
  pressure in hPa and temperature in K.
  '''
  return mr * 100. / satmixwat(p, t)

def mr_to_rh_ice( mr,  p,  t) :
  '''
  Returns relative humidity in %, given the mixing ratio over ice in g/kg,  
  pressure in hPa and temperature in K.
  '''
  return mr * 100. / satmixice(p, t)

def satmix( p, t) :
  '''
  Returns saturation mixing ratio in g/kg, given pressure in hPa and
  temperature in K.
  '''
  twat = t.astype(np.float)
  tice = t.astype(np.float)
  
  twat[twat <= 253.0] = np.nan
  tice[tice >  253.0] = np.nan

  satmixw = satmixwat(p, twat)
  satmixi = satmixice(p, tice)
  
  satmixw[np.isnan(satmixw)] = 0
  satmixi[np.isnan(satmixi)] = 0

  return (satmixw + satmixi)


def satmixwat( p,  t) :
  '''
  Returns saturation mixing ratio over water, in g/kg, given pressure in hPa and
  temperature in K.
  '''
  #es = svpwat(t)
  es = svp(t)
  return (622. * es)/p

def satmixice( p, t) :
  '''
  Returns saturation mixing ratio over ice, in g/kg, given pressure in hPa and
  temperature in K.
  '''
  #es = svpice(t);
  es = svp(t)
  return (622. * es) / p;


def svpwat(t) :
    '''
    Returns saturation vapor pressure over water, in hPa, given temperature in K.
    '''
    a0 =  0.999996876e0
    a1 = -0.9082695004e-2
    a2 =  0.7873616869e-4
    a3 = -0.6111795727e-6
    a4 =  0.4388418740e-8
    a5 = -0.2988388486e-10
    a6 =  0.2187442495e-12
    a7 = -0.1789232111e-14
    a8 =  0.1111201803e-16
    a9 = -0.3099457145e-19
    b = 0.61078e+1
    t -= 273.16
    return (b / ((a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9)))))))))**8.)) 

def svp(t):
    if np.isnan(t):
        return t

    e00 = 611.21
    t00 = 273.16
    ti  = t00 - 23.

    esw = e00 * math.exp(17.502 * (t - t00) / (t-32.19))
    esi = e00 * math.exp(22.587 * (t - t00) / (t+0.7))

    ppv = None
    if t > t00:
        ppv = esw                                          # water phase 
    elif t > ti and t <= t00:
        ppv = esi + (esw - esi)*((t - ti)/(t00 - ti))**2   # mixed phase
    elif t <= ti:
        ppv = esi                                          # ice phase 

    ppv = ppv / 100.0                     # conversion from [pascal] to [mb]

    return ppv


def svpice( t) :
  '''
  Returns saturation vapor pressure over ice, in hPa, given temperature in K.
  The Goff-Gratch equation (Smithsonian Met. Tables,  5th ed., pp. 350, 1984)
  '''
  a = 273.16 / t
  exponent = -9.09718 * (a - 1.) - 3.56654 * log10(a) + 0.876793 * (1. - 1./a) + log10(6.1071)

  return 10.0**exponent
