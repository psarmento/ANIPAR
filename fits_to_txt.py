# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:44:17 2016

@author: psarmento
"""
from astropy.io import fits
from pyfits import Column
import pyfits as pyfs
import numpy as np
from math import *
from shutil import *
import astropy.units as u
import os

wvl0=4.179
wvldelta=6*10**-6

loc = 'fits_files/' #source for spectra, make sure this folder contains ONLY spectra, delete any extra files in the folder
    
stars = os.listdir(loc)

dest = 'raw_spectra/' #destination for the txt format spectra
src = ''

convert = True #if true, the code will convert the vacuum wavelength of APOGEE spectra into air wavelength 
#this is necessary if the user wants to run the next steps of the normalization + synthetic spectra matching

#This simple routine converts spectra from fits files (as they are downloaded from ASPCAP) into a file readable by iSpec

###################################
#Important: if any downloaded spectra has a name with a format different than: spec_apStar-r8-2M00182256+4401222 (+ _ind or _glo).fits.txt,
#change it into that format by including the 2Mass ID of the star in the relevant place
#Otherwise the following routines may FAIL
##################################

#reading the spectra
def read_spec(filename):
    sp = fits.open(filename)
    wavelength0=np.arange(8575, dtype='float')
    wavelength = []
    for w in wavelength0:
        wavelength.append(10**(wvl0+wvldelta*w))
    wavelength = wavelength * u.AA
    flux = sp[1].data
    error = sp[2].data
    return wavelength, flux,error
   
#finding the 2Mass id in the spectra's name
def namefinder(aspcap):
    ap = False
    apogee = ''
    for c in aspcap:
        if c == '2' and ap == False:
            apogee += c
            ap = True
        elif ap == True:
            if c == '.':
                ap = False
            else:
                apogee += c
    return apogee

################################################################################
#air-vacuum conversion

def v2a(wave):
    wave = wave/1000 #converting units into micron
    a = 0
    b1 = 5.792105*10**-2
    b2 = 1.67917*10**-3
    c1 = 238.0185
    c2 = 57.362
    s1 = 1 + a + b1/(c1-1.0/wave**2)+b2/(c2-1.0/wave**2)
    wave2 = wave / s1
    wave2 *= 1000 #converting back into nm
    return wave2

if __name__ == '__main__':

    for n,name in enumerate(stars):
        
        wl, flx,err = read_spec(loc + name)
        
        star_id = namefinder(name)
                
        filenew = 'waveobs flux err\n'
        print(star_id)
        
        if len(flx) in np.arange(0,50):
            for i in range(2):
                for w in range(len(wl)):
                    wave = wl[w]
                    if convert == True:
                        filenew += str(v2a(round(float(wave.base),1)/10)) + ' ' + str(flx[i][w]) + ' ' + str(err[i][w]) + ' \n'
                    else:
                        filenew += str(round(float(wave.base),1)/10) + ' ' + str(flx[i][w]) + ' ' + str(err[i][w]) + ' \n'
                if i == 0:
                    weight = '_ind_weit' #differentiating between individual and global weighting for stars with 2 available spectra
                elif i == 1:
                    weight = '_glo_weit'
                if convert == True:
                    f = open(dest + 'spec_' + name + weight +'_v2a.txt', 'w')
                else:
                    f = open(dest + 'spec_' + name + weight +'.txt', 'w')
                f.write(filenew)
                f.close()
        else:
            for w in range(len(wl)):
                wave = wl[w]
                if convert == True:
                        filenew += str(v2a(round(float(wave.base),1)/10)) + ' ' + str(flx[w]) + ' ' + str(err[w]) + ' \n'
                else:
                        filenew += str(round(float(wave.base),1)/10) + ' ' + str(flx[w]) + ' ' + str(err[w]) + ' \n'
                        
            if convert == True:
                f = open(dest + 'spec_' + name + '_v2a.txt', 'w')
            else:
                f = open(dest + 'spec_' + name + '.txt', 'w')
                
                    
            f.write(filenew)
            f.close()