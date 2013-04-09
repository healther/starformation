# -*- coding: utf-8 -*-
from __future__ import print_function, division
from StringIO import StringIO
import sys
from time import time
import distribution as dist
import starformation
import analysis
import numpy as np
from astropy.io import fits


def main(quiet = False):
    '''main(quiet = False)

This script produces a grid of expected numbers of stars according to the selection 
    criteria of Yusef-Zedah et al. 2009, 702,178-225 The Astrophysical Journal.
    The grid is in av for visual extinction, apera for aperature size and age for the maxage 
    of the starformation size

Parameters
----------
quiet   boolean  if true suppresses all standard output


Returns
----------
A number of fits-files with the sampled stars for different parameters to be specified in
this file.
Standard output is used to report progress, it will print out the parameter set to be 
progressed next and the completeness of the script as
AV aperaturesize maxage completeness ETA
ETA is the time to complete in sec based on the single last operation, huge differences are
to be expected.
'''
    t0 = time()         #timing possibility

    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout


    print(t0,file=output_stream)
    

    sfr = .01
    # star mass function
    def f(x):               # Kouper IMF
    #http://adsabs.harvard.edu/abs/2001MNRAS.322..231K
        if x<.01:
            return 0
        elif x < .08:
            return 62.46192*x**-.3
        elif x < .5:
            return 1.413352*x**-1.8
        elif x < 1.:
            return x**-2.3          # also value 2.7 given eq.(6) or eq (2)
        else:
            return x**-2.3
    f = np.vectorize(f)
    mf = dist.distribution(f, .1, 50.)

    # star formation history
    def g(x):
        return sfr
    g = np.vectorize(g)
    ages = np.logspace(5,7,7)
    sf = [dist.distribution(g, 1000., ages[i]) for i in range(len(ages))]

    t1 = time()                 # finished reading the distributions
    print(t1,file=output_stream)


# setting up model data
    aperas = np.logspace(2, 5, 4)
    avs = np.linspace(10.0, 50.0, 5)
    l = 1
    mpold, tmpnew = 0., time()
    parameters = []
    for i in range(len(avs)):
        for j in range(len(aperas)):
            for k in range(len(ages)):
                tmpold, tmpnew = tmpnew, time()
                starformation.main(massfunction = mf, starformationhistory = sf[k], \
                    A_v = avs[i], sfr = sfr, apera = aperas[j], maxage = ages[k], \
                    appendix = "%s_%03d_%06d_%09d" % ('sim',avs[i],aperas[j],ages[k]), quiet=True)
                print(avs[i],aperas[j],ages[k], l/len(avs)/len(aperas)/len(ages), (len(avs)*len(aperas)*len(ages)-l)*(tmpnew-tmpold),file=output_stream)
                l = l+1
                
                parameters.append([avs[i],aperas[j],ages[k]])

    t2 = time()                 # end of simulation
    print(t2, t1, t2-t1)
    
    print ('number of simulations run: %s' %l , file=output_stream)  
    head = ['#','AV', 'Aperature_size', 'Age']
    f = open('out/__head', 'w')
    f.write( ','.join(head)+'\n' )
    np.savetxt(f, parameters)
    f.close()

    t3 = time()                 # end of saving data

    analysis.main('out')
    print ('analysis complete' , file=output_stream)  
    
    t4 = time()                 # end of analysing data



    print( 'starting script at %f'  %(t0), file=output_stream)
    print( 'initializing       %f'  %(t1-t0), file=output_stream)
    print( "running simulation %f"  %(t2-t1), file=output_stream)
    print( "writing data       %f"  %(t3-t2), file=output_stream)
    print( "analysing data     %f"  %(t4-t3), file=output_stream)
    print( "________________________", file=output_stream)
    print( "total runtime      %f"  %(t4-t0), file=output_stream)
    print( "finishing script   %f"  %t4, file=output_stream)


main()