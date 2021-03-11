# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 18:00:54 2016

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

#This code runs iSpec+Turbospectrum on a sample of normalized stars defined by an input table
#it outputs both the best fit spectrum + table with parameters for characterized stars

##################
#input options

inputs = True

if inputs == True:
    free_params = ["teff", "logg", "MH",'R','vmic','vsini'] #list of free parameters for the code    
    
    L = 'fusion_Mdwarf_New' #Linelist to use, can't actually be changed here due to the way Turbospectrum works
    Code = 'Turbospectrum' #code to use
    #Do not change these lines!
        
    Using_Normalized_stars = True
    
    loggRanger = 0.3 #specifications on the allowed range for each parameter. 
    #Subtract 0.1 from logg and [M/H], 100K from Teff due to boundary conditions inherent to iSpec
    #Actual limits would therefore be from 4.5 to 4.9 when using loggRanger = 0.3 and initial value of 4.7
    teffRanger = 450
    mhRanger = 0.45
    
    #output filename
    filename0 = '7stars_LoggRange_' + str(loggRanger-0.1) + '_TeffRange_' + str(teffRanger-100) + '_MHRange_' + str(mhRanger-0.1)
    
    #Input table name
    Input_Table = 'Mdwarf_BestTemplates.fits'

    loc = ''
    tmp_dir = None
    if Using_Normalized_stars == True:
    
        prefix = ''
     
        folder = prefix + 'normalized_spectra/' #default location for input spectra
        
        save_path = loc + prefix + 'default_output/' #default location for output best fit spectra
    
    norm = False #Set to true if you want to renormalize the spectrum, not recommended
    
    linemasklist = ['NewList_3200'] #linemask to use, can select multiple and code iterates over all of tthem
    
    codes = ['turbospectrum']

################################################################################
#--- iSpec directory -------------------------------------------------------------
ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
#ispec_dir = '/home/virtual/shared/iSpec/'

sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


#--- Change LOG level ----------------------------------------------------------
#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

################################################################################
######################################################
#identifying 2mass id of star used to create spectrum

def starfinder(name):
    out = ''
    for j,c in enumerate(name):
        if c == '2' and out == '' and name[j+1] == 'M':
            out += c
        elif (c != '.' and c != '_') and out != '':
            out += c
        elif (c == '.' or c == '_') and out != '':
            return out
        
################################################################################
#air-vacuum conversion
            
def convert_v2a(star_spectrum):
    star_spectrum_air = vacuum_to_air2(star_spectrum)
    return star_spectrum_air
    
def vacuum_to_air2(spectrum):
    """
    It converts spectrum's wavelengths from vacuum to air
    """
    converted_spectrum = ispec.create_spectrum_structure(spectrum['waveobs'], spectrum['flux'], spectrum['err'])
    wave2 = np.power(spectrum['waveobs'], 2)
    # Compute conversion factor
    a = 0
    b1 = 5.792105*10**-2
    b2 = 1.67917*10**-3
    c1 = 238.0185
    c2 = 57.362
    s1 = 1 + a + b1/(c1-1.0/spectrum['waveobs']**2)+b2/(c2-1.0/spectrum['waveobs']**2)
    converted_spectrum['waveobs'] = spectrum['waveobs'] / s1
    return converted_spectrum

######################################################
# fetching linemask

def create_segments_around_linemasks(star_spectrum, spectral_type = 'Mdwarf'):
    #---Create segments around linemasks -------------------------------------------
    if spectral_type == 'NewList_3200':
        linefile = 'linemask.txt'
    line_regions0 = ispec.read_line_regions(linefile)
    line_regions = ispec.adjust_linemasks(star_spectrum, line_regions0, max_margin=1.0)
    return line_regions
    
#Estimating SNR of studied spectrum
def estimate_snr_from_flux(star_spectrum,points):
    ## WARNING: To compare SNR estimation between different spectra, they should
    ##          be homogeneously sampled (consider a uniform re-sampling)
    #--- Estimate SNR from flux ----------------------------------------------------
    logging.info("Estimating SNR from fluxes...")
    num_points = points
    estimated_snr = ispec.estimate_snr(star_spectrum['flux'], num_points=num_points)
    return estimated_snr

###############################################
#Determining parameters + best fit spectra

def determine_astrophysical_parameters_using_synth_spectra(star_spectrum,line_list,segments,code="synthe",initial_teff = 5050,initial_logg = 4.46,initial_MH = -0.02,initial_vsini = 8, linelist='Classic',normalization= False,weight = '',spt = 'Mdwarf',starname=''):
    #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Radial Velocity determination with template -------------------------------
    logging.info("Radial velocity determination with template...")
    
    if Using_Normalized_stars == False:
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
    #model = "Splines" # "Polynomy"
    #degree = 3
    #nknots = None # Automatic: 1 spline every 5 nm
    from_resolution = to_resolution

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    #median_wave_range=0.05
    #max_wave_range=1.0
    
    if normalization == True:
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.02, model="Fixed value")
        
    #--- Normalize -------------------------------------------------------------
    #normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
                                    
    elif normalization == False:
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        
    normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
    
    #norm_filename = "Normalized_%s_%s_%s.fits" % (starname,weight,m)
    
    #ispec.write_spectrum(normalized_star_spectrum, save_path + norm_filename)

    #--- Model spectra ----------------------------------------------------------
    # Parameters

    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
    initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)

    initial_limb_darkening_coeff = 0.6
    initial_R = to_resolution
    initial_vrad = 0
    max_iterations = 20
    
    initial_alpha = ispec.determine_abundance_enchancements(initial_MH)
    
    logg_range = [initial_logg-loggRanger,initial_logg+loggRanger]
    teff_range = [initial_teff-teffRanger,initial_teff+teffRanger]
    mh_range = [initial_MH-mhRanger,initial_MH+mhRanger]

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "input/atmospheres/MARCS.GES/"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"

    if linelist == 'New':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD+HIT/atomic_lines.tsv"
    elif linelist == 'Classic':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    elif linelist == 'FeH':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/V+HR+FEH/atomic_lines.tsv"
    elif linelist == 'Apogee':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/APOGEE_lines/atomic_lines.tsv"
    elif linelist == 'APOGEE_sergi':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/APOGEE.1500_1700nm/atomic_lines.tsv"
    elif linelist == 'Fusion_list':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Fusion_list/atomic_lines.tsv"
    elif linelist == 'fusion_Mdwarf_New':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Fusion_Mdwarfs_FINAL/atomic_lines.tsv"
    elif linelist == 'souto_lines':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Souto_list/atomic_lines.tsv"
    elif linelist == 'Fusion_FeH_new':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Fusion_FeH_v6/atomic_lines.tsv"
    elif linelist == 'fusion_depth':
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Fusion_depth/atomic_lines.tsv"  
    
    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
    #atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

    element_name = "Fe"
    free_abundances = ispec.create_free_abundances_structure([element_name], chemical_elements, solar_abundances)
    free_abundances['Abund'] += initial_MH # Scale to metallicity

    # Free individual element abundance
    free_abundances = None
    linelist_free_loggf = None
       
    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
            initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
            linemasks=line_list, \
            enhance_abundances=True, \
            use_errors = True, \
            vmic_from_empirical_relation = True, \
            vmac_from_empirical_relation = True, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code,variableRanges = True,logg_range_input=logg_range,teff_range_input=teff_range,mh_range_input=mh_range) #actually running Turbospectrum)
                
                
    t2 = round(params['teff'],1) #output parameter rounding for name of output file
    l2 = round(params['logg'],1)
    m2 = round(params['MH'],1)
    
    ##--- Save results -------------------------------------------------------------
    logging.info("Saving results...")
    #dump_file = "Synth_%s_T_%s_%s_lg_%s_%s_mh_%s_%s_%s_%s_%s_norm=%s_%s.dump" % (code, t, t2, l, l2, m, m2, starname,vs,L,str(Using_Normalized_stars),weight)
    logging.info("Saving results...")
    #ispec.save_results(dump_file, (params, errors, abundances_found, loggf_found, status, stats_linemasks))
    # If we need to restore the results from another script:
    #params, errors, abundances_found, loggf_found, status, stats_linemasks = ispec.restore_results(dump_file)

    logging.info("Saving synthetic spectrum...")
    #Output name for spectra, you can customize it here
    synth_filename = "%s_%s_%s_mask_%s_T_%s_%s_lg_%s_%s_mh_%s_%s_LoggRange%s.fits" % (code,starname,weight,spt,initial_teff, t2, initial_logg, l2, initial_MH, m2,loggRanger)
    ispec.write_spectrum(modeled_synth_spectrum, save_path + synth_filename)
    return params,status,errors



    
####################################################################################
#Creating a fits file to save data

def create_fits():    
    ####   create table
    
    columns = ('APOGEE_ID','T_ASPCAP','log g_ASPCAP','M/H_ASCPAP','vsini_ASPCAP','SNR_ASPCAP','SNR_estimate','Code','Linemask Used','params','T_I','log g_I','MH_I','linelist','T_O','log g_O','MH_O','vmic_O','vmac_O','vsini_O','limb_O','R_O','rms','nsynthesis','wchisq','niter','rchisq','chisq','rwchisq','T_error','logg_error','MH_error')

    dtypes = ('S20','f8','f8','f8','f8','f8','f8','S20','S20','S50', 'f8','f8','f8','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')
    
    MasterTable = Table(names = columns, dtype = dtypes)
    
    return MasterTable

#Opening input table
def open_Stars():
    hdulist = fits.open(Input_Table)
    tbdata = hdulist[1].data
    return tbdata

##################################################
#getting star name + template parameters from name of spectrum file  
    
def extract_Params(star):
    starname = starfinder(star)
    mh = float(MH_finder(star))+0.0
    lg = float(Logg_finder(star))
    tef = float(Teff_finder(star))
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
        
#########################################
# This code creates the output table and oversees the best fit routines
            
def parameters(star):
    starname,t, m, l = extract_Params(star) #getting template parameters from spectrum name
    star_spectrum = ispec.read_spectrum(folder + star) #opening normalized spectra
    
    star_params = [starname,0,0,0,0,0, 'ASPCAP']
    
    #Resample spectrum
    wavelengths = np.arange(np.min(star_spectrum['waveobs']), np.max(star_spectrum['waveobs']), 0.001)
    star_spectrum = ispec.resample_spectrum(star_spectrum, wavelengths, method="linear", zero_edges=True)
    
    #---------------------
    
    line = star_params
    for i in range(25): #creating the new line to add to big results table
        line.append(0)
    line[6] = estimate_snr_from_flux(star_spectrum,5)       
    line[11] = l
    line[9] = str(free_params)
        
    run = True #set to false if you just want to test your setup, set to True to actually calculate Best fits
    print ('Will work on star:')
    print( starname)
    
    
    line[7] = Code
    print ('_____________________________')
    print ('Using code ', Code)
    print( 'Using star ', starname)
    print( '_____________________________')
    
    for spt in linemasklist: #iterating over linemasks
        line[8] = spt
        marg = 1.0
        line_regions = create_segments_around_linemasks(star_spectrum, spectral_type = spt)#,min_dpt,max_dpt, resolution)

        segments = ispec.create_segments_around_lines(line_regions, margin=marg)
        segments.sort()
        
        line[10] = t
        line[12] = m
        line[13] = L
        print( '----------------------------------------------------------------------------')
        print( 'Starting new run')
        print ('Spectral type of linemask used - ', spt)
        print ('Temp - ', t)
        print ('logg - ', l)
        print ('metals - ', m)
        print ('Linelist - ', L)
        #print segments
        #print line_regions
        if run == True:
            params,status,errors = determine_astrophysical_parameters_using_synth_spectra(star_spectrum,line_regions,segments = segments,code=Code,initial_teff = t,initial_logg = l,initial_MH = m,initial_vsini = 0,linelist = L,normalization=norm,spt = spt,starname = starname)
            line[14] = params['teff']
            line[15] = params['logg']
            line[16] = params['MH']
            line[17] = params['vmic']
            line[18] = params['vmac']
            line[19] = params['vsini']
            line[20] = params['limb_darkening_coeff']
            line[21] = params['R']
            line[22] = status['rms']
            line[23] = status['nsynthesis']
            line[24] = status['wchisq']
            line[25] = status['niter']
            line[26] = status['rchisq']
            line[27] = status['chisq']
            line[28] = status['rwchisq']
            line[29] = errors['teff']
            line[30] = errors['logg']
            line[31] = errors['MH']
            #line[34] = cat
        
        MasterTable.add_row(line)
    temp_line_filename = open(save_path + star + '_results_exo.txt', "w")

    temp_line_filename.write(str(line))
    
    temp_line_filename.close()
    return line
    
##############################################################################
    
if __name__ == '__main__': #This is the part that actually runs everything, separated to allow for parallelization of the whole process
    
    print ('Starting by opening Aspcap')
    #Aspcap = open_aspcap(True)
    print ('Creating fits')
    filename = filename0 + '.fits'
    
    MasterTable = create_fits() #creating big table to house results
    
    StarsList = open_Stars() #opening table with spectra that will be characterized
    
    stars = []
    for star in StarsList:
        stars.append('Normalized_' + str(star['APOGEE_ID']) + '_' + str(int(star['T_Norm'])) + '.0_' + str(star['MH_Norm'])+ '_' + str(star['logg_Norm']) + '.fits')
        #this step creates names of spectra that will be characterized as per indications from StarsList about their templates
        #Do not change names of spectra so that this step will run smoothly
    stars.sort()
    
    p = Pool(4) #parallelizing the process, set your number of parallel threads for the code here
    #Set to 1 to run everything in series
    
    table = p.map(parameters, stars)
    for line in table:
        MasterTable.add_row(line)
        
    MasterTable.write(loc + filename,overwrite=True)  #writing big table to file 