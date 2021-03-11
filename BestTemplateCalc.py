#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 10:44:45 2020

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

#This small code finds the best templates for each star, according to the ChiSq values calculated by ChiSqCalc.py
#it can output not only the best template for each star, but also every template within a specified margin of ChiSq value

loc = ''

Output_filename = 'Mdwarf_BestTemplates.fits' #output table
Input_filename = 'Mdwarf_ChiSq_Table.fits' # Input table

margin = 1.10 #margin, set to 1.0 to include only best template per star
#margin will be multiplied by the lowest ChiSq to find all templates within those values

def open_aspcap():
    hdulist = fits.open(Input_filename) #input table
    tbdata = hdulist[1].data
    return tbdata

def create_fits():    
    ####   create table
    columns = ('APOGEE_ID','T_Norm','MH_Norm','logg_Norm','ChiSqMask','ChiSqFull')
    
    dtypes = ('S20','f8','f8','f8','f8','f8')
    
    MasterTable = Table(names = columns, dtype = dtypes)
    
    return MasterTable

BigTable = open_aspcap()
NewTable = create_fits()

BigTable.sort()

CharStars = []
BestChiSq = []
BestFits = []

for line in BigTable: #finding best combination per star
    if line['APOGEE_ID'] not in CharStars:
        BestFits.append(line)
        CharStars.append(line['APOGEE_ID'])
        BestChiSq.append(line['ChiSqMask'])
    elif line['ChiSqMask'] < BestChiSq[-1]:
        BestFits[-1] = line
        CharStars[-1] = line['APOGEE_ID']
        BestChiSq[-1] = line['ChiSqMask']
        
for line in BigTable: #adding other combinations within the margin used
    star = line['APOGEE_ID']
    for k,comp_star in enumerate(CharStars):
        if comp_star == star:
            if line['ChiSqMask'] < BestChiSq[k]*margin and line['ChiSqMask'] != BestChiSq[k]:
                BestFits.append(line)
        
for line in BestFits:
    NewTable.add_row(line)
    
NewTable.write(loc + Output_filename,overwrite=True) #writing results to file