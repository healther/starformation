# -*- coding: utf-8 -*-
from __future__ import print_function, division

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.mlab import griddata
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


import zero


def main(aperture = 1, visual_ex = 0, age=0, bapera=True, bav=True, bage=True):
    """main(aperture = 1, visual_ex = 0, age=0, bapera=True, bav=True, bage=True)

    Parameters
    ----------
    aperture  int:
        defines for which aperture the AV-age-sfr plot is taken
    visual_ex  int:
        defines for which AV the age-aperture-sfr plot is taken
    age        int:
        defines for which age the av-aperture-sfr plot is taken
    bapera     bool:
        if false surpresses the av-age plot
    bav        bool:
        if false surpresses the age-aperture plot
    bage       bool:
        if false surpresses the av-aperture plot

    Returns
    ----------
    A 3d-surface-plot at `the aperture+1`th aperture and a plot at the `visual_ex+1`th AV value
    of the data in 'out/__expected_number'
    """

    avs, aperas, ages, sfr, headers, data = init()

    numbers = np.clip(data[:,3].reshape(len(avs), len(aperas), len(ages)),1.,1000000.) #prevent inf if no stars are detected
    if bapera:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.view_init(45,135)
        X, Y = np.meshgrid(avs, np.log10(ages))
        Z = 559.*sfr/np.transpose(numbers[:,aperture,:])
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
        ax.set_title('Starformation as a function of\n visual extinction Av and age at apperaturesize=%s' % aperas[aperture])

        plt.savefig('plot/3dav-age.svg')

    if bav:
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
        ax2.set_title('Starformation as a function of\napperaturesize and age at Av=%s' % avs[visual_ex])

        plt.savefig('plot/3dage-apera.png')

    if bage:
        fig3 = plt.figure()
        ax3 = fig3.gca(projection='3d')
        ax3.view_init(45,-135)

        X3, Y3 = np.meshgrid(np.log10(aperas), avs)
        Z3 = 559.*sfr/numbers[:,:,age]
        surf3 = ax3.plot_surface(X3, Y3, Z3, rstride=1, cstride=1, cmap=cm.coolwarm)
        cbar3 = fig3.colorbar(surf3)
        cbar3.set_label('starformation rate in M_sun/year')

        x3 = X3.reshape(numbers.shape[0]*numbers.shape[1])
        y3 = Y3.reshape(numbers.shape[0]*numbers.shape[1])
        z3 = Z3.reshape(numbers.shape[0]*numbers.shape[1])   
        ax3.scatter(x3,y3,z3)
        ax3.set_xlabel('log(apera)')
        ax3.set_ylabel('AV')
        ax3.set_zlabel('starformation rate in M_sun/year')
        ax3.set_title('Starformation as a function of\napperaturesize and AV at age=%s years' % ages[age])

        plt.savefig('plot/3dav-apera.png')




def cmd(folder, av, apera, age, color1 = "I4", color2 = "M1", corrected=True, old=False, color='mass'):
    """cmd(folder, av, apera, age, color1 = "I4", color2 = "M1") - creates a CMD

    Parameters
    ----------
    folder   String:
        folder in which the datafile is to be found
    av       float:
        value for the visual extinction as in the filename
    apera    float:
        value of the aperture size as in the filename 
    age      float:
        value of the age as in the filename
    color1   String:
        band filter to be used for the first color 
    color2   String:
        band filter to be used for the second color

    Returns
    ----------
    A color-magnitude-diagram of the data in the 'folder/sim_av_apera_age' fits-file 
    using the `color1-color2` vs `color2`
    """
    xmin = -1.
    xmax = 8.
    ymin = 15.
    ymax = -1.
    selectmax = 0.
    selectmin = 8.
 
    hdulist = fits.open('%s/%s' %(folder,'sim_%03d_%06d_%09d'%(av,apera,age)))
    data = hdulist[1].data
    z = data[color]
    

    if old:
        x = -2.5*(np.log10(data['cflux %s' % color1]/zero.zero_mag[color1]) - np.log10(data['cflux %s' % color2]/zero.zero_mag[color2]))
        y = -2.5*(np.log10(data['cflux %s' % color2]/zero.zero_mag[color2]))
    else:
        if corrected:
            x = -2.5*(np.log10(data['c%s' % color1]/zero.zero_mag[color1]) - np.log10(data['c%s' % color2]/zero.zero_mag[color2]))
            y = -2.5*(np.log10(data['c%s' % color2]/zero.zero_mag[color2]))
        else:
            x = -2.5*(np.log10(data['%s' % color1]/zero.zero_mag[color1]) - np.log10(data['%s' % color2]/zero.zero_mag[color2]))
            y = -2.5*(np.log10(data['%s' % color2]/zero.zero_mag[color2]))
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    scat = ax.scatter(x,y, c=z)
    cb = fig.colorbar(scat)
    cb.set_label('log(%s)' % color)
    ax.plot([1., 4. ], [10., 0.])
    ax.plot([-1., 8.],[13., 4.])
    ax.plot([xmin, xmax],[selectmin, selectmin])
    ax.plot([xmin, xmax],[selectmax, selectmax])

    ax.set_xlabel('%s-%s' % (color1,color2))
    ax.set_ylabel(color2)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title('cmd plot for %s_%s_%s' %(av,apera,age))

    plt.savefig('plot/%s_%s_%s_%s.png' %(av,apera,age,color))


def plot_2d(aperture = 1, visual_ex = 0, age = 0, bapera=True, bav=True, bage=True):
    """plot_2d(aperture = 1, visual_ex = 0, age = 0, bapera=True, bav=True, bage=True) - creates two 2d contour plots

    Parameters
    ----------
    aperture  int:
        defines for which aperture the AV-age-sfr plot is take
    visual_ex  int:
        defines for which AV the age-aperture-sfr plot is take
    age        int:
        defines for which age the av-aperture-sfr plot is taken
    bapera     bool:
        if false surpresses the av-age plot
    bav        bool:
        if false surpresses the age-aperture plot
    bage       bool:
        if false surpresses the av-aperture plot

    Returns
    ----------
    A plot at `the aperture+1`th aperture and a plot at the `visual_ex+1`th AV value
    of the data in 'out/__expected_number'
    """

    avs, aperas, ages, sfr, headers, data = init()

    numbers = np.clip(data[:,3].reshape(len(avs), len(aperas), len(ages)),1.,1000000.) #prevent inf if no stars are detected

    if bapera:
        X, Y = np.meshgrid(avs, np.log10(ages))
        Z = 559.*sfr/np.transpose(numbers[:,aperture,:])
      
        x = X.reshape(numbers.shape[0]*numbers.shape[2])
        y = Y.reshape(numbers.shape[0]*numbers.shape[2])
        z = Z.reshape(numbers.shape[0]*numbers.shape[2])

        fig = plt.figure()
        ax = fig.gca()

        cs = ax.contourf(X,Y,Z)
        #cs1 = ax.contour(X,Y,Z,[0.01,.011,.012,.013,.014,.015],colors=['red','red','red','red','red','red'])
        cs1 = ax.contour(X,Y,Z,levels=[0.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.375,.4], colors=['red','red','red','red','red','red','red','red','red','red','red','red','red'])
        ax.clabel(cs1, fontsize=10, inline=1)
        cbar = plt.colorbar(cs)
    
        ax.scatter(x,y)
        ax.set_xlabel('av')
        ax.set_ylabel('log(age)')
        cbar.set_label('starformation rate in M_sun/year')
        ax.set_title('Starformation as a function of visual extinction Av\nand age at apperaturesize=%s' % aperas[aperture])

        plt.savefig('plot/2dav-age-%s.png' % aperas[aperture])

    if bav:
        fig2 = plt.figure()
        ax2 = fig2.gca()

        X2, Y2 = np.meshgrid(np.log10(ages), np.log10(aperas))
        Z2 = 559.*sfr/numbers[visual_ex,:,:]
        cs2 = ax2.contourf(X2, Y2, Z2)
        cs21 = ax2.contour(X2,Y2,Z2,levels=[0.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.375,.4], colors=['red','red','red','red','red','red','red','red','red','red','red','red','red'])
        ax2.clabel(cs21, fontsize=10, inline=1)
        cbar2 = plt.colorbar(cs2)

        x2 = np.log10(data[:,2].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
        y2 = np.log10(data[:,1].reshape(len(avs), len(aperas), len(ages))[visual_ex,:,:].reshape(numbers.shape[1]*numbers.shape[2]))
        ax2.scatter(x2,y2)
        ax2.set_xlabel('log(age)')
        ax2.set_ylabel('log(apera)')
        cbar2.set_label('starformation rate in M_sun/year')
        ax2.set_title('Starformation as a function of\n apperaturesize and age at Av=%s' % avs[visual_ex])

        plt.savefig('plot/2dage-apera-%s.png' % avs[visual_ex])

    if bage:
        fig3 = plt.figure()
        ax3 = fig3.gca()

        X3, Y3 = np.meshgrid(np.log10(aperas), avs)
        Z3 = 559.*sfr/numbers[:,:,age]
        cs3 = ax3.contourf(X3, Y3, Z3)
        cs31 = ax3.contour(X3, Y3, Z3,levels=[0.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.375,.4], colors=['red','red','red','red','red','red','red','red','red','red','red','red','red'])
        ax3.clabel(cs31, fontsize=10, inline=1)

        cbar3 = fig3.colorbar(cs3)
        cbar3.set_label('starformation rate in M_sun/year')

        x3 = X3.reshape(numbers.shape[0]*numbers.shape[1])
        y3 = Y3.reshape(numbers.shape[0]*numbers.shape[1])
        ax3.scatter(x3,y3)
        ax3.set_xlabel('log(apera)')
        ax3.set_ylabel('AV')
        ax3.set_title('Starformation as a function of\napperaturesize and AV at age=%s years' % ages[age])
        cbar3.set_label('starformation rate in M_sun/year')

        plt.savefig('plot/2dav-apera-%s.png' % ages[age])

def init():
    """init() - aquires the base data

    """
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


def histogram(folder, av, apera, age):
    """histogram(folder, av, apera, age, color1 = "I4", color2 = "M1")

    Parameters
    ----------
    folder   String:
        folder in which the datafile is to be found
    av       float:
        value for the visual extinction as in the filename
    apera    float:
        value of the aperture size as in the filename 
    age      float:
        value of the age as in the filename

    Returns
    ----------
    A histogram of the massdistribution in logspace with 20 bins from -.5 to 1.5 M_sun
    of all the sampled stars and restricted to the selected stars 
    """
    
    hdulist = fits.open('%s/%s' %(folder,'sim_%03d_%06d_%09d'%(av,apera,age)))
    data = hdulist[1].data

    fig = plt.figure()
    ax = fig.add_subplot(111)


    pdf, bins, patches = ax.hist(data['mass'], range=(-.5,1.5), label='all', bins=20)
    dat = data['mass'][data['sel'] == 1]
    pdf, bins, patches = ax.hist(dat, range=(-.5,1.5), label='selected', bins=20)
    ax.legend()
    #ax.set_yscale('log')
    ax.set_ylabel('number of stars')
    ax.set_xlabel('log(mass/M_sun)')
    ax.set_title('mass distribution simulated\nwith AV=%s, age=%s years, apertursize=%s AU' % (av,age,apera))

    plt.savefig('plot/hist-%s-%s-%s.png' % (av,apera,age))










