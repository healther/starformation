# -*- coding: utf-8 -*-
from __future__ import print_function, division
from StringIO import StringIO
from time import time
import sys

import numpy as np
from astropy.io import fits

import distribution as dist
import starformation
import analysis
import functions


def main(n = 150000, quiet = False):
    """main(n = 150000, quiet = False)

    This script produces a grid of expected numbers of stars according to the selection 
        criteria of Yusef-Zedah et al. 2009, 702,178-225 The Astrophysical Journal.
        The grid is in av for visual extinction, apera for aperature size and age for the maxage 
        of the starformation size

    Parameters
    ----------
    n       integer:
        number of stars to be sampled per parameter set
    quiet   boolean:
        if true suppresses all standard output


    Returns
    ----------
    A number of fits-files with the sampled stars for different parameters to be specified in
    this file.
    Standard output is used to report progress, it will print out the parameter set to be 
    progressed next and the completeness of the script as
    AV aperaturesize maxage completeness ETA
    ETA is the time to complete in seconds based on the single last operation
    """
    t0 = time()         #timing possibility

    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout


    print(t0,file=output_stream)
    

    sfr = .01
    # star mass function
    kroupa = np.vectorize(functions.kroupa)
    mf = dist.Distribution(kroupa, .1, 50.)

    #star formation history
    constant_sfr = np.vectorize(functions.constant_sfr)
    
    ages = np.logspace(5,7,7)
    sf = [dist.Distribution(constant_sfr, 1000., ages[i]) for i in range(len(ages))]
    #sfr = [150000*mf.mean()/(ages[i]-1000.) for i in range(len(ages))]

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
                    A_v = avs[i], sfr = n, apera = aperas[j], maxage = ages[k], \
                    appendix = "%s_%03d_%06d_%09d" % ('sim',avs[i],aperas[j],ages[k]), quiet=True, precise=False)
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