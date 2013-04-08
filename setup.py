# -*- coding: utf-8 -*-
from __future__ import print_function
import urllib2
import os
from StringIO import StringIO
import sys
import tarfile



def main(quiet=False):
    '''setup.main(quiet=False)

Pulling all necessary files for using the starformation script

This script pulls the fits-files for the radiation models from the MPIA-folder of Thomas
Robitaille, if the files have been moved feel free to contact robitaille@mpia.de

It also pulls the extinction law from 
http://caravan.astro.wisc.edu/protostars/info.php?topic=sedfitter_results 
if you have problems with them please contact
http://caravan.astro.wisc.edu/protostars/contact.php
'''
    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout

    newpath = r'%s/models' % os.getcwdu()
    if not os.path.exists(newpath): os.makedirs(newpath)
    newpath = r'%s/out' % os.getcwdu()
    if not os.path.exists(newpath): os.makedirs(newpath)
    existing = sorted(os.listdir('%s/%s' % (os.getcwdu(), 'models'))) 

    urls = [
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/parameters.fits.gz',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/2J.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/2H.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/2K.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/I1.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/I2.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/I3.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/I4.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/M1.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/M2.fits',
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/M3.fits',
            'http://caravan.astro.wisc.edu/protostars/files/extinction_law.tar.gz'
          ]
    file_names = [
            'models/parameters.fits.gz',
            'models/2J.fits',
            'models/2H.fits',
            'models/2K.fits',
            'models/I1.fits',
            'models/I2.fits',
            'models/I3.fits',
            'models/I4.fits',
            'models/M1.fits',
            'models/M2.fits',
            'models/M3.fits',
            'models/extinction_law.tar.gz']

    for i in range(len(urls)):
        if not os.path.isfile(file_names[i]):
            f = open(file_names[i], 'wb')
            f.write(urllib2.urlopen(urls[i]).read())
            f.close()
            print('Downloaded %s from %s' % (file_names[i],urls[i]), file=output_stream)

    if not os.path.isfile('modesl/extinction_law.ascii')
        f = tarfile.open('models/extinction_law.tar.gz', 'r:gz')
        try: f.extractall()
        finally: f.close()

