# -*- coding: utf-8 -*-
from __future__ import print_function
from time import time
from StringIO import StringIO
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
from astropy.table import Table, Column
from astropy.io import fits

import distribution as dist  
import functions


def main(massfunction = 0, starformationhistory = 0, A_v = 10.0, sfr = .01, apera = 24000,\
 maxage = 2000000., distance = 8.0, appendix='default', quiet=0, precise=0):
    """main(massfunction = 0, starformationhistory = 0, A_v = 10.0, sfr = .01, apera = 24000,\
          maxage = 2000000., distance = 8.0, appendix='default', quiet=0, precise=0)

    Creates a sample of stars

    Parameters
    ----------
    massfunction            distribution:
        relatively in mass, with lower and upper restriction, see also what the distribution must provide
    starformation history   distribution:
        relatively in age, with lower and upper restriction, see also what the distribution must provide
    A_v       float:
        value for the visual extinction 
    sfr       float:
        average star formation rate in M_sun/year (only if precise = True)
    apera     float:
        used aperture size for selecting the fluxes of the protostars
    maxage    float:
        age of the star formation site, sfr is assumed to be constant
    distance  float:
        distance to the simulated starformation site
    appendix  String:
        sets the outputfilename, default is the starting time (via time.time())
    quiet     boolean:
        if true (=1) suppresses all standard output
    precise   boolean:
        if true (=1) sample single star till expected mass reached based on the 
        integrated starformationhistory times the starformationrate
        else sfr is the number of expected stars and the starformationrate is
        calculated by the cumulated mass divided by the formation time

    The distributions must provide an object which has the following members:
        float    cdf(float x)   returns the integrated distribution up to x, is used to calculate
                                the expected mass
        float    _upperbound    returns the upper limit of the distribution, is used to calculate
                                the expected mass
        float[]  sample(int n)  returns an array of n floats, sampled from the distribution
        float    mean()         returns the mean value of the distribution


    Returns
    ----------
    returns a fits file in the out-folder, either using the appendix as filename or the time of the
          starting of the script in order to prevent overwriting existing files
          In the header of the fits-file are the values: A_v, sfr, apera, maxage and distance recorded
          In the data part are the age, mass, modelnumber and the uncorrected and corrected fluxes
    """
    
    if quiet:
        output_stream = StringIO()
    else:
        output_stream = sys.stdout

    t0 = time()                 
    if appendix=='default':  # making sure not to overwrite former output
        appendix=t0          # by using the starting time as an unique id
    #parameter settings
    k_v = 211.4   # opacity in v_band in cm^2/g
    # wavelength of the corresponding filterband in microns
    wavelength = [1.235, 1.662, 2.159, 3.550, 4.493, 5.731, 7.872, 23.68, 71.42, 155.9] 
    models = ['2H', '2J', '2K', 'I1', 'I2', 'I3', 'I4', 'M1', 'M2', 'M3']


    if massfunction == 0 and starformationhistory == 0:
        # star mass function
        kroupa = np.vectorize(functions.kroupa)
        massfunction = dist.Distribution(kroupa, .1, 50.)

        #star formation history
        constant_sfr = np.vectorize(functions.constant_sfr)
        starformationhistory = dist.Distribution(constant_sfr, 1000., maxage)


    cumass = 0.  #sampled mass
    stars = []  #storing for the sample
    sfh = starformationhistory

    t1 = time()  #startup completed

    if precise:
        n = 0
        exmass = sfh.cdf()(sfh._upperbound)*sfr     #expected mass formed
        while cumass < exmass:
            mass, age = massfunction.sample(), sfh.sample()
            cumass = cumass + mass
            stars.append([n, age, mass])
            if n % 10000 == 0:
                print (n, cumass, file=output_stream)                                 #reporting progress
            n = n+1
    else:
        n = sfr
        mass, age = massfunction.sample(n), sfh.sample(n)
        cumass = np.sum(mass)
        exmass = n * massfunction.mean()
        stars = [[i, age[i], mass[i]] for i in range(n)]
    sfr = cumass/(sfh._upperbound-sfh._lowerbound)  #average star formation rate

    print ('number of sampled stars: %s' %n , file=output_stream)  
    print ('mass of sampled stars: %s' % cumass , file=output_stream)  
    print ('mean mass: %s' % (cumass/n), file=output_stream)
    print ('expected mass of stars: %s' % exmass , file=output_stream)
    t2 = time()  # sampleing completed


    # python code for model contact
    #initial parameters
    model = [ fits.open('models/%s.fits' % mod) for mod in models ]    # fits-data for the model
    param = fits.open('models/parameters.fits.gz')  # modelparameter
    app_num = [ np.interp(apera, model[i][2].data.field(0), range(model[i][2].data.field(0).size)) for i in range(len(models)) ] 


    # sampling viewing angle
    angle = np.random.random_integers(0,9,len(stars))
    #reading model grid
    mass = param[1].data['MASSC'][::10]
    age = param[1].data['TIME'][::10]
    grid = np.vstack([age, mass]).transpose()

    #converting to logspace
    stars = np.asarray(stars)
    grid = np.log10(grid)
    stars[:,1:] = np.log10(stars[:,1:])

    output = stars.tolist()  #creating output
    
    #normalizing for nearest neighbor search
    grid[0,:] = grid[0,:]/(grid[0,:].max() - grid[0,:].min())
    grid[1,:] = grid[1,:]/(grid[1,:].max() - grid[1,:].min())
    stars[1,:] = stars[1,:]/(grid[0,:].max() - grid[0,:].min())
    stars[2,:] = stars[2,:]/(grid[1,:].max() - grid[1,:].min())

    t3 = time()  #model data load complete

    tree = scipy.spatial.cKDTree(grid,leafsize=10)  #search tree
    matches = [tree.query(star[1:] , k=1)[1] for star in stars]  #saves matches with (dist, index)

    t4 = time()  #matching sample to data complete

    # extracting fluxes
    fluxes = [0 for j in range(len(models)) ]
    indices = 10*np.asarray(matches) + angle
    for j in range(len(models)):
        fluxes[j] = model[j][1].data[indices]['TOTAL_FLUX'][:,app_num[j]]



    # applying extinction
    extinction = np.loadtxt('models/extinction_law.ascii')
    k_lambda = np.interp(wavelength, extinction[:,0], extinction[:,1])
    correctionfactor = 10.**(-.4 * A_v * k_lambda / k_v)

    newfluxes = [0 for j in range(len(models)) ]
    for j in range(len(models)):
        newfluxes[j] = np.asarray(fluxes[j]) * correctionfactor[j] * (1./distance)**2


    t5 = time()  #extracting fluxes complete

    # saving data
    fluxes = np.asarray(fluxes)
    newfluxes = np.asarray(newfluxes)
    output = np.vstack([np.asarray(output).transpose(), matches, fluxes, newfluxes]).transpose()

    # create table
    # data table
    t = Table()
    t.add_column(Column(name='age', data=output[:,1]))
    t.add_column(Column(name='mass', data=output[:,2]))
    t.add_column(Column(name='model', data=output[:,3]))
    for i in range(len(models)):
        t.add_column(Column(name='%s' % models[i], data=output[:,4+i]))
    for i in range(len(models)):
        t.add_column(Column(name='c%s' % models[i], data=output[:,4+len(models)+i]))
    # head table
    header = Table()
    header.add_column(Column(name='AV', data = [A_v]))
    header.add_column(Column(name='SFR', data = [sfr]))
    header.add_column(Column(name='APPERA', data = [apera])  )   
    header.add_column(Column(name='MAXAGE', data = [maxage]))
    header.add_column(Column(name='DIST', data = [distance]))


    fits.writeto('out/%s' % appendix, np.array(t), clobber=True)
    fits.append('out/%s' % appendix, np.array(header), clobber=True)
    
    t6 = time()  #saving complete

    # timing possibility for optimization efforts

    print( 'starting script at %f'  %(t0), file=output_stream)
    print( 'initializing       %f'  %(t1-t0), file=output_stream)
    print( "sampleing          %f"  %(t2-t1), file=output_stream)
    print( "model data load    %f"  %(t3-t2), file=output_stream)
    print( "matching model     %f"  %(t4-t3), file=output_stream)
    print( "extracting fluxes  %f"  %(t5-t4), file=output_stream)
    print( "saving             %f"  %(t6-t5), file=output_stream)
    print( "________________________", file=output_stream)
    print( "total runtime      %f"  %(t6-t0), file=output_stream)
    print( "finishing script   %f"  %t6, file=output_stream)

#main(sfr = .08)  # for testing purposes and directly called from bash