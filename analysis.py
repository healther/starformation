# -*- coding: utf-8 -*-
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os
from decimal import Decimal
from StringIO import StringIO
import sys
from astropy.io import fits
from astropy.table import Table, Column



def main(folder, quiet=0):
    '''analysis.main(folder, quiet=0)

This script gives the number of "observed" stars from the sampled datafiles in "folder"/ 
according to the selection criteria from Yusef-Zadeh et al

Parameters
----------
folder  String   Specifiecs the folder where the files are
quiet   boolean  =1 suppresses all standard output

Returns
-------
returns a file named __expected_number in which contains a numpy array of the simulation parameters
number, A_v, Aperature_size, Age and the expected 'detected' number
The first line of the file is an ','-seperated head of the contained informations
'''

    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout



    color1 = "I4"             #filter system for first color of CMD
    color2 = "M1"             #filter system for second color of CMD
    min_mag = 8.              #minimal observation limit
    max_mag = 0.              #maximal observation limit

#getting file list
    files = sorted(os.listdir('%s/%s' % (os.getcwdu(), folder))) 
    out = []

    for fil in files:
#only using files created by the automated simulation
        if fil.startswith('sim_') and not 'settings' in fil.encode("ascii"):
            print ("%s/%s" % (folder,fil.encode("ascii")), file=output_stream)
      

      # Read in
            hdulist = fits.open('%s/%s' %(folder,fil))
            data = hdulist[1].data

      #calculating magnitudes from fluxes and converting to CMD-data
            x = -2.5*(np.log10(data['c%s' % color1]/64130) - np.log10(data['c%s' % color2]/7140))
            y = -2.5*(np.log10(data['c%s' % color2]/7140))

            
            sel = np.logical_and( (y > -10./3. * (x-1.) + 10.), np.logical_and(max_mag < y, y < min_mag))
            sel = np.logical_and(sel, y < -x + 12.)
            n = sum(sel)
            t = Table(hdulist[1].data)
            if 'sel' in t.columns:
                t.remove_column('sel')
            t.add_column(Column(name='sel', data=sel.astype('int')))
            
            hdulist[1].data = np.array(t)
            tmp, av, apera, age = fil.split('_')
            #hdulist = fits.open('%s/%s' %(folder,fil), 'write')
            #hdulist.writeto('%s/%s' %(folder,fil), clobber=True)
            fits.update('%s/%s' %(folder,fil), np.array(t), ext = 1, clobber=True)
            out.append([Decimal(av), Decimal(apera), Decimal(age), n])

    #writing obtained data to "folder/__expected_number"
    head = ['#', 'AV', 'Aperature_size', 'Age', 'Expected_number']
    f = open('%s/__expected_number' % folder, 'w')
    f.write(','.join(head)+'\n' )
    np.savetxt(f, out)
    f.close()
   
    print ("Analysed %s files and saved output to %s" % (len(out),'%s/__expected_number' % folder), file=output_stream)
