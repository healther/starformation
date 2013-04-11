#! /usr/bin/env python
# -*- coding: utf-8 -*-

import scipy.integrate as inte
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpol



class Distribution( object ):
    """distribution(func, a, b)

    A distribution object represents a probability distribution, the provided function is taken as the
    probability density function and an cumulative density function and its inverse are calculated. `a` and
    `b` are taken as the lower and the upper bound respectively.
    `func` need not to be normalized

    Parameters
    ----------
    func    function:
        specifies the probability density function of the distribution
    a       float:
        specifies the lower bound of the distribution
    b       float:
        specifies the upper bound of the distribution

    Members
    ----------
    pdf():          
        returns the probability density function, takes a float as an argument
    cdf():         
        returns the cumulative probability density function, takes a float as an argument
    ppf():          
        returns the inverse cumulative probability funciton, takes a float as an argument
    norm():
        returns the normalisation factor (integral over pdf)
    sample(int n):
        returns an array of n sampled values
    mean():
        returns the mean value of the distribution
    _lowerbound:
        lowerbound of the distribution
    _upperbound:
        upperbound of the distribution
    """
    def __init__(self, func = lambda x: 1, a = 0., b = 1000.):
        self._pdf = func
        self._lowerbound = a
        self._upperbound = b
	self.cdf_init(self._lowerbound, self._upperbound)
	self._norm = self.cdf()(self._upperbound)
        t, step = np.linspace(self._lowerbound, self._upperbound, 10000, endpoint=True, retstep=True)
        y = [step*( self.cdf()(t[i]) + self.cdf()(t[i+1]) )/2. for i in (range(t.size - 1))]
        self._mean = self._upperbound - sum(y)/self._norm

    def pdf(self):
        return self._pdf

    def cdf_init(self, a, b):
        t, step = np.linspace(a, b, 10000, endpoint=True, retstep=True)
        y = np.array(t)
        y[0], i = 0, 0
        for i in (range(t.size - 1)):
            y[i+1] = y[i] + step*( self.pdf()(t[i]) + self.pdf()(t[i+1]) )/2
        self._x, self._y = t, y
    
    def cdf(self):
        return interpol.interp1d(self._x, self._y)

    def ppf(self):                              # inverse cdf 
        return interpol.interp1d(self._y, self._x)

    def norm(self):
        return self._norm

    def sample(self, n):
        t = np.random.uniform(0,self.norm(),n)
        return self.ppf()(t)
    
    def sample(self):
        t = np.random.uniform(0,self.norm())
        return self.ppf()(t)

    def mean(self):                     
        return self._mean



