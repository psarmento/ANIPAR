#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:51:21 2020

@author: psarmento
"""

import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
import ispec

#This is a simple code to calculate ChiSq similarities between normalized spectra and the synthetic spectra used to normalize them
#Run normalizer_Mdwarfs_grid.py beforehand in order to get the normalized spectra themselves

loc = ''

filename = 'Mdwarf_ChiSq_Table.fits' #output table's file name

linemaskname = 'linemask.txt'

folder = 'normalized_spectra/' #source for normalized spectra

stars = os.listdir(loc + folder)

def extract_Params(star):
    starname = starfinder(star)
    mh = str(float(MH_finder(star))+0.0)
    lg = Logg_finder(star)
    tef = Teff_finder(star)
    return starname, tef, mh, lg
    
def Teff_finder(star):
    under_counter = 0
    t = ''
    for c in star:
        if under_counter == 2 and len(t)<4:
            t += c
        if c == '_':
            under_counter += 1
    return t

def Logg_finder(star):
    under_counter = 0
    point_counter = 0
    t = ''
    for c in star:
        if c == '_':
            under_counter += 1
        elif under_counter == 4 and point_counter<3:
            if c == '.':
                point_counter += 1
            if point_counter < 2:
                t += c
    return t

def MH_finder(star):
    under_counter = 0
    point_counter = 0
    t = ''
    for c in star:
        if c == '_':
            under_counter += 1
        elif under_counter == 3 and point_counter<3:
            if c == '.':
                point_counter += 1
            if point_counter < 2:
                t += c
    return t

def starfinder(name):
    out = ''
    for j,c in enumerate(name):
        if c == '2' and out == '' and name[j+1] == 'M':
            out += c
        elif (c != '.' and c != '_') and out != '':
            out += c
        elif (c == '.' or c == '_') and out != '':
            return out

def open_template(tef, mh, lg):
    match0 = 'Grid/grid_synth_'
    match1 = match0 + str(tef) + '.0_' + str(lg) + '_' + str(mh) + '.fits'
    synth_spectrum = ispec.read_spectrum(match1)
    return synth_spectrum

def resample_spectrum(star_spec):
    #--- Resampling  --------------------------------------------------------------
    logging.info("Resampling...")
    wavelengths = np.arange(1514, 1694, 0.02)
    resampled_star_spectrum = ispec.resample_spectrum(star_spec, wavelengths, method="bessel", zero_edges=True)
    return resampled_star_spectrum

def create_comparing_mask(waveobs, linemasks, segments=None):
    # Build wavelength points from regions
    wfilter = None
    for region in linemasks:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        # Consider only lines that are inside segments
        if segments is not None:
            in_segment1 = np.logical_and(segments['wave_base'] <= wave_base, segments['wave_top'] >= wave_base)
            in_segment2 = np.logical_and(segments['wave_base'] <= wave_top, segments['wave_top'] >= wave_top)
            in_segment = np.logical_and(in_segment1, in_segment2)
            if np.all(in_segment == False):
                continue

        if wfilter is None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_linemask = np.zeros(len(waveobs))
    waveobs_linemask[wfilter] = True # Consider fluxes only for selected line masks

    return waveobs_linemask

def ChiSqCalc(sp1,sp2,mask=None):
    if mask == None:
        residuals = sp1['flux'] - sp2['flux']
    else:
        residuals = sp1['flux'][mask] - sp2['flux'][mask]
    ChiSq = np.sum((residuals)**2)
    return ChiSq

def create_fits():    
    ####   create table
    
    columns = ('APOGEE_ID','T_Norm','MH_Norm','logg_Norm','ChiSqMask','ChiSqFull')
    
    dtypes = ('S20','f8','f8','f8','f8','f8')
    
    MasterTable = Table(names = columns, dtype = dtypes)
    
    return MasterTable

def create_segments_around_linemasks(star_spectrum):
    #weight = weightfinder(starname)
    #---Create segments around linemasks -------------------------------------------
    linefile = linemaskname        
    line_regions0 = ispec.read_line_regions(linefile)
    line_regions = ispec.adjust_linemasks(star_spectrum, line_regions0, max_margin=1.0)
    return line_regions

BigTable = create_fits()

run = True #set this to false to test the functions themselves, without running the full code

if run == True:
    
    for i,star in enumerate(stars): #running across all normalized spectra
        
        print( star, i)
        
        starname, tef, mh, lg = extract_Params(star) #getting the template parameters from the name of the normalized spectra
        template = open_template(tef, mh, lg) #finding the template used to normalize the spectra
        spectra = ispec.read_spectrum(loc + folder + star) #opening the template
        
        linemasks = create_segments_around_linemasks(spectra) #applying the linemask to the spectra
        
        template = resample_spectrum(template) #resampling both template and normalized spectra for comparison
        spectra = resample_spectrum(spectra)
        
        waveobs = spectra['waveobs']
        comparing_mask = create_comparing_mask(waveobs, linemasks) #masking the compared regions according to linemask
        mask = np.where(comparing_mask==1.0)
        
        ChiSqMask = ChiSqCalc(spectra,template,mask) #calculating Chi Sq in linemask regions
        ChiSqFull = ChiSqCalc(spectra,template) #calculating Chi Sq in every wavelength region
        
        line = [starname, tef, mh, lg, ChiSqMask, ChiSqFull] #creating a line summarizing the results
        BigTable.add_row(line) #adding it to table
    
BigTable.write(loc + filename,overwrite=True) #writing table to file
    
    