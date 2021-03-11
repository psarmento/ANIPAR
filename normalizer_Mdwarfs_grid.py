# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 18:00:54 2016

@author: psarmento
"""

import os
import sys
import numpy as np
import logging
from multiprocessing import Pool
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
import scipy.signal
import string

#this code normalizes all spectra in selected folder

##################
#input options

inputs = True

if inputs == True:
    
    loc = ''#location for the normalized spectra and output, can be customized, using own folders by default
    
    filename = 'Normalized_Mdwarfs_data.fits' #data for normalized M dwarfs will be stored here
    
    #txt_folder = 'raw_txt/'
    txt_folder = 'raw_spectra/' #location of raw spectra to normalize
    
    #stars = os.listdir(loc + 'norm_stars_data/')
    stars = os.listdir(loc + txt_folder)
    
    norm = "Template" #renormalize the spectra?

################################################################################
#--- iSpec directory -------------------------------------------------------------
ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
#ispec_dir = '/home/virtual/shared/iSpec/'
#save_path = loc + 'default_output/' 
save_path = loc + 'normalized_spectra/'
#location where output spectra will be saved

sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

#--- Change LOG level ----------------------------------------------------------
#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

###################################################################################
#Getting the isochrone values

iso_file = np.loadtxt('isochrones.dat')

teffs = [[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
loggs = [[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]]]

age = iso_file[0][2]
age0 = iso_file[0][2]
ages = [age]

metal = iso_file[0][1]
metal0 = iso_file[0][1]
metals = [metal]

c = 0
r = 0

Combos = []

for k in iso_file:
    Teff = (10**k[7])/1000
    logg = k[8]
    if 3.0 < Teff < 4.000 and logg>4.5:
        NewCombo = [k[1]+0.0,1000*round(Teff,1),round(logg,1)]
        if NewCombo not in Combos:
            Combos.append(NewCombo)


######################################################
#identifying combination + 2mass id of star used to create spectrum
    
def weightfinder(star):
    if 'glo' in star:
        return 'glo'
    elif 'ind' in star:
        return 'ind'
    
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
            elif c in string.digits or c == '+' or c == '-' or c == 'M':
                apogee += c
    return apogee
    
######################################################
#fetching linemask
        
def create_segments_around_linemasks(star_spectrum):
    #weight = weightfinder(starname)
    #---Create segments around linemasks -------------------------------------------
    linefile = 'linemask.txt'    #this is the linemask used, you can customize it here    
    line_regions0 = ispec.read_line_regions(linefile)
    line_regions = ispec.adjust_linemasks(star_spectrum, line_regions0, max_margin=1.0)
    return line_regions

################################################
#actual normalization code, based on ispec, using its functionalities
    
def normalise_spectra(star_spectrum,star_params,normalization= False,weight = '', cat = 0, v2a = False,starname=''):
    #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Radial Velocity determination with template -------------------------------
    logging.info("Radial velocity determination with template...")
    if v2a == True:
        star_spectrum = convert_v2a(star_spectrum)
    #rv = np.round(models[0].mu(), 2) # km/s
    #rv_err = np.round(models[0].emu(), 2) # km/s
    #--- Radial Velocity correction ------------------------------------------------
    #logging.info("Radial velocity correction... %.2f +/- %.2f" % (rv, rv_err))
    #star_spectrum = ispec.read_spectrum('/home/psarmento/Documents/iSpec_treated_spectra/txt/txt_stars3/spec_aspcapStar-r5-v603-'+starname+'.fits.txt')
    #star_spectrum = ispec.correct_velocity(star_spectrum, rv)
    #--- Resolution degradation ----------------------------------------------------
    # NOTE: The line selection was built based on a solar spectrum with R ~ 47,000 and VALD atomic linelist.
    #from_resolution = 80000
    to_resolution = 22000
    #star_spectrum = ispec.convolve_spectrum(star_spectrum, to_resolution, from_resolution)
    #--- Continuum fit -------------------------------------------------------------
    model = "Polynomy" # "Splines" 
    degree = 2
  #  nknots = None # Automatic: 1 spline every 5 nm
    from_resolution = to_resolution

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
  #  median_wave_range=0.05
  #  max_wave_range=1.0
    
    if normalization == True:
        star_continuum_model = ispec.fit_continuum(star_spectrum, independent_regions=segments,
                                                   from_resolution=from_resolution, \
                                    nknots=nknots, degree=degree, \
                                    median_wave_range=median_wave_range, \
                                    max_wave_range=max_wave_range, \
                                    model=model, order=order, \
                                    automatic_strong_line_detection=True, \
                                    strong_line_probability=0.5, \
                                    use_errors_for_fitting=True)
                                    
    elif normalization == False:
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        
    if normalization == 'Template':
        best_fit,matches = match_spectrum(star_params)
        synth_spectrum = ispec.read_spectrum(best_fit)
    
        #--- Continuum fit -------------------------------------------------------------
        model = "Template"
        nknots = None # Automatic: 1 spline every 5 nm (in this case, used to apply a gaussian filter)
        from_resolution = 22500
        median_wave_range=5.0
    
        #strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
        #strong_lines = "/home/psarmento/Documents/iSpec_treated_spectra/Solar_apogee/New/NEWLINES/linemasks_22500/solar_apogee_fixed.dat_star_fitted_linemasks_turbospectrum_Fusion_list.txt"
        strong_lines = None
        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                    ignore=strong_lines, \
                                    nknots=nknots, \
                                    median_wave_range=median_wave_range, \
                                    model=model, \
                                    template=synth_spectrum)
    print (star_spectrum['err'])
    print (matches)
                                    
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    
    teffy = matches[0]
    mety = star_params[3]
    logy = star_params[2]
    
    norm_filename = "Normalized_%s_%s_%s_%s.fits" % (starname,teffy,mety,logy) #name for the newly normalized spectrum
    #DO NOT change the names of output spectra, as they are needed for consistency reasons in other codes
    
    clean_spectra = clean_spectrum(normalized_star_spectrum)
    
    print (clean_spectra['err'])
    
    #new_snr = estimate_snr_from_flux(clean_spectra,5)
    
    #final_spectra = accordion_spectra(clean_spectra)
    final_spectra = clean_spectra
    ispec.write_spectrum(final_spectra, save_path + norm_filename)
    return matches,final_spectra, synth_spectrum

##############################################
#fetching template spectrum
    
def match_spectrum(star_params,fixed=False):
    teff = float(star_params[1])
    logg = star_params[2]
    mh = star_params[3]
    match0 = 'Grid/grid_synth_'

   # if mann == 0:
   #     photo_teff = visual
   # else:
   #     photo_teff = (visual + mann)/2
    print (teff, mh, logg)#, visual,mann, photo_teff
    
    if fixed == False:
        teff_match = str(int(teff)) + '.0'
        mh_match = str(mh)
        logg_match = str(logg)
        match1 = match0 + teff_match + '_' + logg_match + '_' + mh_match + '.fits'
        matches = [teff_match,mh_match,logg_match,teff]
    else:
        teff_match = '3200'
        mh_match = '0.4'
        logg_match = '5.4'
        match1 = match0 + teff_match + '_' + logg_match + '_' + mh_match + '.fits'
        matches = [teff_match,mh_match,logg_match,photo_teff]
        
    return match1,matches
#############################################
#cleaning poor parts of normalized spectrum + estimating SNR
    
def clean_spectrum(star_spectrum):
    #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Clean fluxes and errors ---------------------------------------------------
    logging.info("Cleaning fluxes and errors...")
    flux_base = 0.0 #any flux values below this will be fixed at selected value
    flux_top = 1.02 #any flux values above this will be fixed at selected value
    err_base = 0.0
    err_top = 1.0
    #star_spectrum = spectrum_filter(star_spectrum)
    ffilter = (star_spectrum['flux'] > flux_base)
    #efilter = (star_spectrum['err'] > err_base) & (star_spectrum['err'] <= err_top)
    #wfilter = np.logical_and(ffilter, efilter)
    clean_star_spectrum = star_spectrum[ffilter]
    topfilter = (clean_star_spectrum['flux'] >= flux_top)
    #print topfilter
    #print clean_star_spectrum[topfilter]
    for i in range(0,len(clean_star_spectrum['flux'])):
        if topfilter[i] == True:
            clean_star_spectrum[i]['flux'] = flux_top
   # print clean_star_spectrum[topfilter]
    #norm_filename = "Clean_%s.fits" % (starname)
    #ispec.write_spectrum(clean_star_spectrum, save_path_clean + norm_filename)
    return clean_star_spectrum
    
def estimate_snr_from_flux(star_spectrum,points):
    ## WARNING: To compare SNR estimation between different spectra, they should
    ##          be homogeneously sampled (consider a uniform re-sampling)
    #--- Estimate SNR from flux ----------------------------------------------------
    logging.info("Estimating SNR from fluxes...")
    num_points = points
    estimated_snr = ispec.estimate_snr(star_spectrum['flux'], num_points=num_points)
    return estimated_snr

####################################################################################
#Creating a fits file to save data

def create_fits():    
    ####   create table
    
    columns = ('APOGEE_ID','T_O','log_g_O','M/H_O','Vis_Teff','Mann_Teff','ASPCAP?','teff_M','MH_M','log_g_M','photo_teff','ChiSq')
    
    dtypes = ('S20','f8','f8','f8','f8','f8','S20','f8','f8','f8','f8','f8')
    
    MasterTable = Table(names = columns, dtype = dtypes)
    
    return MasterTable

##############################################################################
def Normalizing_Stars(star):

    weight = weightfinder(star) #Identification of combination used for combination of spectrum (individual or global)
    run = True
    if weight != 'ind': #can be changed to another combination (glo)
        star_spectrum = ispec.read_spectrum(loc + txt_folder + star) #fetching spectra
        #star_spectrum = convert_v2a(star_spectrum) #converting to air wavelength
        starname = namefinder(star) #retrieving 2mass id from the name of the spectrum
        linemasks = create_segments_around_linemasks(star_spectrum) #matching linemask to stellar spectrum
        
        star_paramys = [star,0,0,0,0,0, 'ASPCAP']
        
        cat = 1
        
        print ('Will work on star:')
        print( star, weight)
        print( starname)
        if run == True:
            for combo in Combos: #running through different combos for spectrum normalization
                star_paramys[1] = combo[1]
                star_paramys[2] = combo[2]
                star_paramys[3] = combo[0]
                print ('Will normalize with template: ',star_paramys)
                matches,spectra, template = normalise_spectra(star_spectrum,star_paramys,normalization=norm, weight = weight, cat = cat, v2a = False,starname=starname)
                print( matches)

#    return line
    
##############################################################################
    
if __name__ == '__main__':
    
    print ('Starting by opening Aspcap')
    #Aspcap = open_aspcap(True)
    print( 'Creating fits')
    
    MasterTable = create_fits()
    
    stars.sort()
    
    p = Pool(10) #this step parallelizes the process, set your number of parallel threads for the code here
    table = p.map(Normalizing_Stars, stars)
            
