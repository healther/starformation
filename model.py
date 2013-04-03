# -*- coding: utf-8 -*-

import numpy as np
import scipy.spatial
from astropy.io import fits

# python code for model contact


#initial parameters
model = fits.open('models/I3.fits')             # fits-data for the model
param = fits.open('models/parameters.fits.gz')
#stars =                                         # star sample
apera = 24000.                                  # aperature size

#


app_num = np.interp(apera, model[2].data.field(0), range(model[2].data.field(0).size))
#if abs(1-(num % 1)) < .1:
    #num = num + 1 - (num % 1)           # for good fits use next aperature size
#else:
#    'to be implemented'                 # for bad fits interpolate

#mass = param[1].data.field('massc')
#age  = param[1].data.field('time')
#name = param[1].data.field('model_name')

grid = [param[1].data[10.*i][1:3] for i in range(param[1].data.size /10.) ] #model data 0: name, 1: age, 2: mass ; only one per model

tree = scipy.spatial.cKDTree(grid,leafsize=10)                   #search tree
matches = [tree.query( star , k=1) for star in stars]           #safes matches with (dist, index)

# extracting fluxes
fluxes = [(model[1].data[match[1]][1][num] , model[1].data[match[1]][2][num]) for match in matches]


