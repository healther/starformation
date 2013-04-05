# -*- coding: utf-8 -*-
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os
from decimal import Decimal
from StringIO import StringIO
import sys



def main(folder, quiet=0):
    '''This script gives the number of "observed" stars from the sampled datafiles in "folder"/ 
according to the selection criteria from Yusef-Zadeh et al

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
            #f = open("%s/%s" % (folder,fil.encode("ascii")), 'r')
            #headers = f.readline().strip().split(',')
            #data = np.loadtxt(f)
            #f.close()
      

      # Read in
            hdulist = fits.open('file.fits')
            av = hdulist[1].header['age']
            data = hdulist[1].data

      ##selecting the relevant columns
            #c1 = headers.index('corrected_flux %s' % color1)
            #c2 = headers.index('corrected_flux %s' % color2)

      #calculating magnitudes from fluxes and converting to CMD-data
            x = -2.5*(np.log10(data['corrected_flux %s' % color1]/64130) - np.log10(data[:,c2]/7140))
            y = -2.5*(np.log10(data['corrected_flux %s' % color2]/7140))


          # efficiency? accuracy?
            n=0
      #selecting "observed" stars
            for i in range(len(x)):
                if (y[i] < -10./3. * (x[i]-1.) + 10.) and (max_mag < y[i] < min_mag):
                    n = n+1
            tmp, av, apera, age = fil.split('_')
            out.append([Decimal(av), Decimal(apera), Decimal(age), n])

    #writing obtained data to "folder/__expected_number"
    head = ['#', 'AV', 'Aperature_size', 'Age', 'Expected_number']
    f = open('%s/__expected_number' % folder, 'w')
    f.write(','.join(head)+'\n' )
    np.savetxt(f, out)
    f.close()
   
    print ("Analysed %s files and saved output to %s" % (len(out),'%s/__expected_number' % folder), file=output_stream)
