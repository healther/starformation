# -*- coding: utf-8 -*-
from __future__ import print_function
from StringIO import StringIO
import sys
import distribution as dist
import starformation
import analysis
import numpy as np
from astropy.io import fits


def main(quiet = False):
    ''' This script produces a grid of expected numbers of stars according to the selection 
    criteria of Yusef-Zedah et al. 2009, 702,178-225 The Astrophysical Journal.
    The grid is in av for visual extinction, apera for aperature size and age for the maxage 
    of the starformation size


'''
    models = ['2H', '2J', '2K', 'I1', 'I2', 'I3', 'I4', 'M1', 'M2', 'M3']
    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout

    sfr = .001

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
    ages = np.linspace(500000, 2000000, 11)
    sf = [dist.distribution(g, 10000., ages[i]) for i in range(10,11)]

# setting up model data
    aperas = np.linspace(15000, 50000, 11)
    avs = np.linspace(5.0, 50.0, 10)
    k = 0.
    parameters = []
    for i in range(len(avs)):
        for i in range(len(aperas)):
            for j in range(len(ages)):
                starformation.main(mf, sf[i], avs[i], .001, app_num[j], ages[i], "%s_%s_%s" % (av,apera,age), quiet=True)
                print(av, apera, age, k/len(avs)/len(aperas)/float(len(ages)), file=output_stream)
                k = k+1
                parameters.append([av,apera,age])

    starformation.main(massfunction = mf, starformationhistory = sf[0], A_v = 10.0, sfr = .001, apera = aperas[2], maxage = 2000000.)
    print ('number of simulations run: %s' %k , file=output_stream)  
    head = ['AV', 'Aperature_size', 'Age']
    f = open('out/__head', 'w')
    f.write( ','.join(head)+'\n' )
    np.savetxt(f, parameters)
    f.close()
    analysis.main('out', quiet=True)
    print ('analysis complete' , file=output_stream)  
 
