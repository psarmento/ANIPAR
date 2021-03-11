#!/usr/bin/env python
#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool

import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
import astropy.units as u

##################

save_loc = 'default_output/' #put your destination for the spectral grid here

#printing to file

import sys

###################################################################################
#Getting the isochrone values

iso_file = np.loadtxt('isochrones.dat') #isochrone sources, can be customizable. remember to delete file header first

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

maxTeff = 4.00 #maximum value for the Teff used in the grid (/1000)
maxLogg = 4.5 #maximum value for the logg used in the grid

Combos = [] #this will contain the combinations used for the grid

#getting isochrone values, rounding to nearest multiples of 100K for Teff, 0.1 dex for [M/H] and logg
for k in iso_file:
    Teff = (10**k[7])/1000
    logg = k[8]
    if Teff < maxTeff and logg>maxLogg:
        NewCombo = [k[1],1000*round(Teff,1),round(logg,1)]
        if NewCombo not in Combos:
            Combos.append(NewCombo)


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

def synthesize_spectrum(teff = 3200,logg = 5.00, MH = 0.0,code="Turbospectrum",name='2mass'):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    #teff = 3200 
    #logg = 5.00
    #MH = 0.00
    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
    macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
    alpha = ispec.determine_abundance_enchancements(MH)
    vsini = 1.60 # Sun
    limb_darkening_coeff = 0.6
    resolution = 22500
    wave_step = 0.01

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    regions = None
    wave_base = 1500 # Magnesium triplet region
    wave_top = 1700


    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "input/atmospheres/MARCS.GES/"#"modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"
    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/Fusion_Mdwarfs_FINAL/atomic_lines.tsv"
    
    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
    #atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])
    #fixed_abundances = ispec.read_solar_abundances(solar_abundances_file) # No fixed abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH}):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print(msg)

    # Enhance alpha elements + CNO abundances following MARCS standard composition
  #  alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = ispec.determine_abundance_enchancements(MH)
  #  alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = (0,0.02,0.0,-0.04)
  #  abundances =  ispec.enhance_solar_abundances(solar_abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

    fixed_abundances = None
    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH}, code=code)

    # Synthesis
    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    print('going to start synthesis with parameters Teff %s MH %s logg %s' % (teff,MH,logg))
    # synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
    #         atmosphere_layers, teff, logg, MH, atomic_linelist, isotopes, solar_abundances, \
    #         fixed_abundances, microturbulence_vel = microturbulence_vel, \
    #         macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
    #         R=resolution, regions=regions, verbose=1,
    #         code=code)
        
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
            fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1,
            code=code)

    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    synth_filename = "grid_synth_%s_%s_%s.fits" % (teff,logg,MH) #default name for each synthetic spectrum, can be customized
    ispec.write_spectrum(synth_spectrum, save_loc + synth_filename)

#for combo in Combos: #running through every parameter combination, creating a synthetic spectrum for each of them
#    synthesize_spectrum(teff = combo[1],logg = combo[2], MH = combo[0])
    
synthesize_spectrum(teff = 3800,logg = 4.6, MH = 0.2)