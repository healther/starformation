# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os
from decimal import Decimal


def main(folder):

    color1 = "I4"             #filter system for first color of CMD
    color2 = "M1"             #filter system for second color of CMD
    min_mag = 8.
    max_mag = 0.

    files = sorted(os.listdir('%s/%s' % (os.getcwdu(), folder)))
    out = []

    for fil in files:
        if not ('settings' in fil.encode("ascii") or fil.startswith('__')):
            f = open("%s/%s" % (folder,fil.encode("ascii")), 'r')
            headers = f.readline().split(',')
            data = np.loadtxt(f)
            f.close()


            c1 = headers.index('corrected_flux %s' % color1)
            c2 = headers.index('corrected_flux %s' % color2)


            x = -2.5*(np.log10(data[:,c1]/64130) - np.log10(data[:,c2]/7140))
            y = -2.5*(np.log10(data[:,c2]/7140))


          # efficiency?
            n=0
            for i in range(len(x)):
                if (y[i] < -3./8. * x[i] + 83./8.) and (max_mag < y[i] < min_mag):
                    n = n+1
            av, apera, age = fil.split('_')
            out.append([Decimal(av), Decimal(apera), Decimal(age), n])


    head = ['#', 'AV', 'Aperature_size', 'Age', 'Expected_number']
    f = open('%s/__expected_number' % folder, 'w')
    f.write(','.join(head)+'\n' )
    np.savetxt(f, out)
    f.close()
   

