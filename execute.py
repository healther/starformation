# -*- coding: utf-8 -*-


import starformation
import analysis
import numpy as np


def main():

        #main(A_v = 10.0, sfr = .001, apera = 24000, age = 2000000., appendix='default'):
    k = 0.
    parameters = []
    for av in np.linspace(5.0, 50.0, 10):
        for apera in np.linspace(15000, 50000, 8):
            for age in np.linspace(500000, 2000000, 12):
                starformation.main(av, .001, apera, age, "%s_%s_%s" % (av,apera,age))
                print av, apera, age, k/10./8./12.
                k = k+1
                parameters.append([av,apera,age])
    head = ['AV', 'Aperature_size', 'Age']
    f = open('out/__head', 'w')
    f.write( ','.join(head)+'\n' )
    np.savetxt(f, parameters)
    f.close()

 
