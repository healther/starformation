# -*- coding: utf-8 -*-
from time import time
import distribution as dist  
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
from astropy.io import fits

def main(A_v = 10.0, sfr = .001, apera = 24000, maxage = 2000000., appendix='default'):
    '''Creates a sample of stars

input:
A_v     float   value for the visual extinction 
sfr     float   Star formation rate in M_sun/year, is assumed to be constant over maxage
apera   float   used aperature size for selecting the fluxes of the protostars
maxage  float   age of the star formation site, sfr is assumed to be constant

output:
returns two files in the folder 'out/' the _settings file contains the used values of 
          A_v, sfr, apera, maxage, number of sampled stars, their cumulated mass and
          the expected mass
'''
    t0 = time()                 
    if appendix=='default':             # making sure not to overwrite former output
        appendix=t0                     # by using the starting time as an unique id
    #parameter settings
    k_v = 211.4             # opacity in v_band in cm^2/g
            # wavelength of the corresponding filterband in microns
    wavelength = [1.235, 1.662, 2.159, 3.550, 4.493, 5.731, 7.872, 23.68, 71.42, 155.9] 
    models = ['2H', '2J', '2K', 'I1', 'I2', 'I3', 'I4', 'M1', 'M2', 'M3']



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
    sf = dist.distribution(g, 10000., maxage)

    # todo - change sample an expected number of stars (expmass/averagemass)
    cumass = 0.                                                     #sampled mass
    exmass = sf.cdf()(sf._upperbound)-sf.cdf()(sf._lowerbound)      #expected mass formed
    stars = []                                               #storing for the sample
    n = 0

    t1 = time()                     # startup completed


    while cumass < exmass:
        mass, age = mf.sample(1)[0], sf.sample(1)[0]
        cumass = cumass + mass
        stars.append([n, age, mass])
#        if n % 10000 == 0:
#            print n, cumass                                 #reporting progress
        n = n+1


    t2 = time()                      # sampleing completed


    # python code for model contact
    #initial parameters
    model = [ fits.open('models/%s.fits' % mod) for mod in models ]    # fits-data for the model
    param = fits.open('models/parameters.fits.gz')  # modelparameter
    app_num = [ np.interp(apera, model[i][2].data.field(0), range(model[i][2].data.field(0).size)) for i in range(len(models)) ] 

    # to do:
        # check for interpolation of aperature size

    # sampling viewing angle
    angle = [np.random.random_integers(9) for i in range(len(stars))]
    #reading model grid
    grid = [param[1].data[10*i][1:3] for i in range(param[1].data.size /10) ] #model data 0: name, 1: age, 2: mass ; only one per model

    #converting to logspace
    stars = np.asarray(stars)
    grid = np.log10(grid)
    stars[:,1:] = np.log10(stars[:,1:])

    output = stars.tolist()                                 #creating output

    #normalizing for nearest neighbor search
    grid[0,:] = grid[0,:]/(grid[0,:].max() - grid[0,:].min())
    grid[1,:] = grid[1,:]/(grid[1,:].max() - grid[1,:].min())
    stars[1,:] = stars[1,:]/(grid[0,:].max() - grid[0,:].min())
    stars[2,:] = stars[2,:]/(grid[1,:].max() - grid[1,:].min())

    t3 = time()                       #model data load complete

    tree = scipy.spatial.cKDTree(grid,leafsize=10)                   #search tree
    matches = [tree.query(star[1:] , k=1)[1] for star in stars]      #saves matches with (dist, index)

    t4 = time()                       #matching sample to data complete

    # extracting fluxes
    fluxes = [0 for j in range(len(models)) ]
    for j in range(len(models)):
    #    fluxes[j] = [ [model[j][1].data[10*matches[i] + angle[i]][1][app_num[j]], model[j][1].data[10*matches[i] + angle[i]][2][app_num[j]] ] for i in range(len(matches)) ]
        fluxes[j] = [ model[j][1].data[10*matches[i] + angle[i]][1][app_num[j]] for i in range(len(matches)) ]


    # applying extinction
    extinction = np.loadtxt('models/extinction_law.ascii')
    k_lambda = np.interp(wavelength, extinction[:,0], extinction[:,1])
    correctionfactor = 10.**(-.4 * A_v * k_lambda / k_v)

    newfluxes = [0 for j in range(len(models)) ]
    for j in range(len(models)):
        newfluxes[j] = np.asarray(fluxes[j]) * correctionfactor[j] 


    t5 = time()                       #extracting fluxes complete

    # saving data to output: #, log10(age), log10(mass), modelmatch, (flux, flux_error, corrected_flux, corrected_flux_error) for each model
    for i in range(len(output)):
        output[i].append(matches[i])
        for j in range(len(models)):
            output[i].append(fluxes[j][i])
            output[i].append(newfluxes[j][i])
    #        output[i].append(fluxes[j][i][0])          possible to use if fluxerrors 
    #        output[i].append(fluxes[j][i][1])          are necessary
    #        output[i].append(newfluxes[j][i][0])
    #        output[i].append(newfluxes[j][i][1])

    # creating the output file
    head = ['#', 'age', 'mass', 'model']
    for mod in models:
        head.append('flux %s' % mod)
        head.append('corrected_flux %s' % mod)
    f = open('out/%s' % appendix, 'w')
    f.write( ','.join(head)+'\n' )
    np.savetxt(f, output)
    f.close()

    # creating the settings file
    f = open('out/%s_settings' % appendix, 'w')
    settings =            '%s #visual extinction    A_v  \n' % A_v
    settings = settings + '%s #star formation rate  sfr  \n' % sfr
    settings = settings + '%s #star formation time  time \n' % maxage
    settings = settings + '%s #aperature size       apera\n' % apera
    settings = settings + '%s #number of sampled stars   \n' % len(stars)
    settings = settings + '%s #cumulated sampled mass    \n' % cumass
    settings = settings + '%s #expected mass             \n' % exmass
    f.write(settings)
    f.close()
    
    t6 = time()                       #saving complete

# timing possibility for optimization efforts

#    print 'starting script at %f'  %(t0)
#    print 'initializing       %f'  %(t1-t0)
#    print "sampleing          %f"  %(t2-t1)
#    print "model data load    %f"  %(t3-t2)
#    print "matching model     %f"  %(t4-t3)
#    print "extracting fluxes  %f"  %(t5-t4)
#    print "saving             %f"  %(t6-t5)
#    print "________________________"
#    print "total runtime      %f"  %(t6-t0)
#    print "finishing script   %f"  %t6

