# -*- coding: utf-8 -*-
from __future__ import print_function, division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def main():
    aperature = 1
    visual_ex = 0


    f = open('out/__expected_number', 'r')
    headers = f.readline().strip().split(',')[1:]
    data = np.loadtxt(f)
    f.close()


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    avs = np.linspace(10.0, 50.0, 5)
    aperas = np.logspace(2, 5, 4)
    ages = np.linspace(500000, 2000000, 11)
    numbers = np.clip(data[:,3].reshape(5, 4, 11),5.,100000.) #prevent inf if no stars are detected


    X, Y = np.meshgrid(avs, ages)
    Z = 559.*.08/np.transpose(numbers[:,aperature,:])
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
    fig.colorbar(surf)


    x = X.reshape(55)
    y = Y.reshape(55)
    z = Z.reshape(55)


    ax.scatter(x,y,z)
    ax.set_xlabel('av')
    ax.set_ylabel('age')
    ax.set_zlabel('starformation rate in M_sun/year')
    #ax.set_zlim(0.,1.)
    ax.set_title('Starformation as a function of\n visual extinction Av and age at apperaturesize=1kpc')

    #plt.savefig('plot/av-age.svg')

    fig2 = plt.figure()
    ax2 = fig2.gca(projection='3d')


    X2, Y2 = np.meshgrid(ages, np.log10(aperas))
    Z2 = 559.*.08/numbers[visual_ex,:,:]
    surf2 = ax2.plot_surface(X2, Y2, Z2, rstride=1, cstride=1, cmap=cm.coolwarm)
    fig2.colorbar(surf2)


    x2 = data[:,2].reshape(5, 4, 11)[visual_ex,:,:].reshape(44)
    y2 = np.log10(data[:,1].reshape(5, 4, 11)[visual_ex,:,:].reshape(44))
    z2 = 559.*.08/data[:,3].reshape(5, 4, 11)[visual_ex,:,:].reshape(44)
    ax2.scatter(x2,y2,z2)
    ax2.set_xlabel('age')
    ax2.set_ylabel('log(apera)')
    ax2.set_zlabel('starformation rate in M_sun/year')
   # ax2.set_zlim(0.,1.)
    ax2.set_title('Starformation as a function of\n visual extinction apperaturesize and age at Av=10')



   # plt.savefig('plot/age-apera.png')




def cmd(folder, av, apera, age):
    color1 = "I4"             #filter system for first color of CMD
    color2 = "M1"             #filter system for second color of CMD
    xmin = -1.
    xmax = 8.
    ymin = 15.
    ymax = -1.
    selectmax = 0.
    selectmin = 8.

    #f = open("%s/sim_%s_%s_%s" % (folder,av,apera,age), 'r')
    #headers = f.readline().strip().split(',')
    #data = np.loadtxt(f)
    #f.close()
  
    hdulist = fits.open('%s/%s' %(folder,'sim_%03d_%06d_%07d'%(av,apera,age)))
#            av = hdulist[1].header['age']
    data = hdulist[1].data
  
    x = -2.5*(np.log10(data['cflux %s' % color1]/64130) - np.log10(data['cflux %s' % color2]/7140))
    y = -2.5*(np.log10(data['cflux %s' % color2]/7140))


    #c1 = headers.index('corrected_flux %s' % color1)
    #c2 = headers.index('corrected_flux %s' % color2)

    #x = -2.5*(np.log10(data[:,c1]/64130) - np.log10(data[:,c2]/7140))
    #y = -2.5*(np.log10(data[:,c2]/7140))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x,y)
    ax.plot([1., 4. ], [10., 0.])
    ax.plot([xmin, xmax],[selectmin, selectmin])
    ax.plot([xmin, xmax],[selectmax, selectmax])


    ax.set_xlabel('%s-%s' % (color1,color2))
    ax.set_ylabel(color2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title('%s_%s_%s.svg' %(av,apera,age))

    plt.savefig('plot/%s_%s_%s.svg' %(av,apera,age))
