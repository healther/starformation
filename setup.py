# -*- coding: utf-8 -*-
from __future__ import print_function
import urllib2
import os




def main(quiet=False):
    '''Pulling all necessary files for using the starformation script

'''
    print(quiet)
    if quiet:
        def print(*args):
            pass

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
            'http://www.mpia-hd.mpg.de/~robitaille/share/andreas/M3.fits'
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
            'models/M3.fits']

    for i in range(len(urls)):
        print('blub')
        if not os.path.isfile(file_names[i]):
            u = urllib2.urlopen(urls[i])
            f = open(file_names[i], 'wb')
            meta = u.info()
            file_size = int(meta.getheaders("Content-Length")[0])
            print ("Downloading: %s Bytes: %s" % (file_names[i], file_size))

            file_size_dl = 0
            block_sz = 8192
            while True:
                buffer = u.read(block_sz)
                if not buffer:
                    break

                file_size_dl += len(buffer)
                f.write(buffer)
                status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
                status = status + chr(8)*(len(status)+1)
                print (status)

            f.close()

    
    