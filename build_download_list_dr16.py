# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:33:16 2016

@author: psarmento
"""
from math import *
import numpy as np

from astropy.io import fits

hdulist2 = fits.open('7stars.fits') #file name containing stars for download
#must be in ASPCAP table format - table with all stars available at https://www.sdss.org/dr14/irspec/spectro_data/

tbdata2 = hdulist2[1].data
x2 = tbdata2['FILE']
y2 = tbdata2['LOCATION_ID']
z2 = tbdata2['TARGET_ID']
t2 = tbdata2['TELESCOPE']
stri = 'https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/' #DR14, must be changed if wanting to download spectra from other DR
ns = 100

for j,name in enumerate(x2):
    name = name[10:]
    loc = y2[j]
    tar = z2[j]
    tel = t2[j]
    if loc != 1:
        stri += 'apo25m/' + str(loc) + '/apStar-r8-' + str(name)
    
    elif 'Mdwarfs' in tar:
        stri += 'apo1m/Mdwarfs/apStar-r8-' + str(name)
    else:
        stri += 'apo1m/calibration/apStar-r8-' + str(name)
        
    
    stri += '\n' + 'https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/stars/'

f = open('down_list_example.txt', 'w') #saving to file
f.write(stri)
f.close()

#after running this code and creating the file, the spectra can be downloaded by opening a command line and executing the code:
#wget -i down_list_example.txt
#The spectra will appear in the folder where the code is executed