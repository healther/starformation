# -*- coding: utf-8 -*-
from __future__ import print_function, division
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def main(aperature = 1, visual_ex = 0):
'''main(aperature = 1, visual_ex = 0)

Parameters
----------
aperature  int    defines for which aperature the AV-age-sfr plot is take
visual_ex  int    defines for which AV the age-aperature-sfr plot is take

Returns
----------
A 3d-surface-plot at `the aperature+1`th aperature and a plot at the `visual_ex+1`th AV value
of the data in 'out/__expected_number'
'''

    #f = open('out/__head', 'r')
    #headers = f.readline().strip().split(',')[1:]
    #data = np.loadtxt(f)
    #f.close()
    #avs = np.unique(data[:,headers.index('AV')])
    #aperas = np.unique(data[:,headers.index('Aperature_size')])
    #ages = np.unique(data[:,headers.index('Age')])
    
    #hdulist = fits.open('%s/%s' %('out','sim_%03d_%06d_%09d'%(avs[0],aperas[0],ages[0])))
    #sfr = hdulist[2].data['SFR'][0]

    #f = open('out/__expected_number', 'r')
    #headers = f.readline().strip().split(',')[1:]
    #data = np.loadtxt(f)
    #f.close()
    avs, aperas, ages, sfr, headers, data = init()

    numbers = np.clip(data[:,3].reshape(len(avs), len(aperas), len(ages)),1.,1000000.) #prevent inf if no stars are detected

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(45,135)
    X, Y = np.meshgrid(avs, np.log10(ages))
    Z = 559.*sfr/np.transpose(numbers[:,aperature,:])
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
    cbar = fig.colorbar(surf)
    cbar.set_label('starformation rate in M_sun/year')

    x = X.reshape(numbers.shape[0]*numbers.shape[2])
    y = Y.reshape(numbers.shape[0]*numbers.shape[2])
    z = Z.reshape(numbers.shape[0]*numbers.shape[2])

    ax.scatter(x,y,z)
    ax.set_xlabel('av')
    ax.set_ylabel('log(age)')
    ax.set_zlabel('starformation rate in M_sun/year')
    ax.set_title('Starformation as a function of\n visual extinction Av and age at apperaturesize=1kpc')

    plt.savefig('plot/3dav-age.svg')

    fig2 = plt.figure()
    ax2 = fig2.gca(projection='3d')
    ax2.view_init(45,45)

    X2, Y2 = np.meshgrid(np.log10(ages), np.log10(aperas))
    Z2 = 559.*sfr/numbers[visual_ex,:,:]
    surf2 = ax2.plot_surface(X2, Y2, Z2, rstride=1, cstride=1, cmap=cm.coolwarm)
    cbar2 = fig2.colorbar(surf2)
    cbar2.set_label('starformation rate in M_sun/year')

    x2 = np.log10(data[:,2].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
    y2 = np.log10(data[:,1].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
    z2 = 559.*sfr/data[:,3].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2])
    ax2.scatter(x2,y2,z2)
    ax2.set_xlabel('log(age)')
    ax2.set_ylabel('log(apera)')
    ax2.set_zlabel('starformation rate in M_sun/year')
    ax2.set_title('Starformation as a function of\napperaturesize and age at Av=10')

    plt.savefig('plot/3dage-apera.png')




def cmd(folder, av, apera, age, color1 = "I4", color2 = "M1"):
'''cmd(folder, av, apera, age, color1 = "I4", color2 = "M1") - creates a CMD

Parameters
----------
folder   String
av       float
apera    float
age      float
color1   String
color2   String

Returns
----------
A color-magnitude-diagram of the data in the 'folder/sim_av_apera_age' fits-file 
using the `color1-color2` vs `color2`
'''
    xmin = -1.
    xmax = 8.
    ymin = 15.
    ymax = -1.
    selectmax = 0.
    selectmin = 8.
 
    hdulist = fits.open('%s/%s' %(folder,'sim_%03d_%06d_%09d'%(av,apera,age)))
    data = hdulist[1].data
  
    x = -2.5*(np.log10(data['c%s' % color1]/64130) - np.log10(data['c%s' % color2]/7140))
    y = -2.5*(np.log10(data['c%s' % color2]/7140))

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

    plt.savefig('plot/%s_%s_%s.png' %(av,apera,age))


def plot_2d(aperature = 1, visual_ex = 0):
'''plot_2d(aperature = 1, visual_ex = 0) - creates two 2d contour plots

Parameters
----------
aperature  int    defines for which aperature the AV-age-sfr plot is take
visual_ex  int    defines for which AV the age-aperature-sfr plot is take

Returns
----------
A plot at `the aperature+1`th aperature and a plot at the `visual_ex+1`th AV value
of the data in 'out/__expected_number'
'''

    #f = open('out/__head', 'r')
    #headers = f.readline().strip().split(',')[1:]
    #data = np.loadtxt(f)
    #f.close()
    #avs = np.unique(data[:,headers.index('AV')])
    #aperas = np.unique(data[:,headers.index('Aperature_size')])
    #ages = np.unique(data[:,headers.index('Age')])
    
    #hdulist = fits.open('%s/%s' %('out','sim_%03d_%06d_%09d'%(avs[0],aperas[0],ages[0])))
    #sfr = hdulist[2].data['SFR'][0]

    #f = open('out/__expected_number', 'r')
    #headers = f.readline().strip().split(',')[1:]
    #data = np.loadtxt(f)
    #f.close()

    avs, aperas, ages, sfr, headers, data = init()

    numbers = np.clip(data[:,3].reshape(len(avs), len(aperas), len(ages)),1.,1000000.) #prevent inf if no stars are detected

    X, Y = np.meshgrid(avs, np.log10(ages))
    Z = 559.*sfr/np.transpose(numbers[:,aperature,:])
   
    x = X.reshape(numbers.shape[0]*numbers.shape[2])
    y = Y.reshape(numbers.shape[0]*numbers.shape[2])
    z = Z.reshape(numbers.shape[0]*numbers.shape[2])

    fig = plt.figure()
    ax = fig.gca()

    cs = ax.contourf(X,Y,Z)
    #cs1 = ax.contour(X,Y,Z,[0.01,.011,.012,.013,.014,.015],colors=['red','red','red','red','red','red'])
    cs1 = ax.contour(X,Y,Z,levels=[0.01,.011,.012,.013,.014,.015],cmap=plt.cm.jet)
    cbar = plt.colorbar(cs)
 
    ax.scatter(x,y)
    ax.set_xlabel('av')
    ax.set_ylabel('log(age)')
    cbar.set_label('starformation rate in M_sun/year')
    ax.set_title('Starformation as a function of visual extinction Av\nand age at apperaturesize=1kpc')

    plt.savefig('plot/2dav-age.png')

    fig2 = plt.figure()
    ax2 = fig2.gca()

    X2, Y2 = np.meshgrid(np.log10(ages), np.log10(aperas))
    Z2 = 559.*sfr/numbers[visual_ex,:,:]
    cs2 = ax2.contourf(X2, Y2, Z2)
    cs21 = ax2.contour(X2,Y2,Z2,[0.01,.011,.012,.013,.014,.015])
    cbar2 = plt.colorbar(cs2)

    x2 = np.log10(data[:,2].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
    y2 = np.log10(data[:,1].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
    ax2.scatter(x2,y2)
    ax2.set_xlabel('log(age)')
    ax2.set_ylabel('log(apera)')
    cbar2.set_label('starformation rate in M_sun/year')
    ax2.set_title('Starformation as a function of\n apperaturesize and age at Av=10')

    plt.savefig('plot/2dage-apera.png')


def init():
'''init() - aquires the base data

'''
    f = open('out/__head', 'r')
    headers = f.readline().strip().split(',')[1:]
    data = np.loadtxt(f)
    f.close()
    avs = np.unique(data[:,headers.index('AV')])
    aperas = np.unique(data[:,headers.index('Aperature_size')])
    ages = np.unique(data[:,headers.index('Age')])
    
    hdulist = fits.open('%s/%s' %('out','sim_%03d_%06d_%09d'%(avs[0],aperas[0],ages[0])))
    sfr = hdulist[2].data['SFR'][0]

    f = open('out/__expected_number', 'r')
    headers = f.readline().strip().split(',')[1:]
    data = np.loadtxt(f)
    f.close()

    return avs, aperas, ages, sfr, headers, data

