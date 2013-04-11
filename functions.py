# -*- coding: utf-8 -*-
 
def kroupa(x):               # Kouper IMF
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

    

    
def constant_sfr(x):
    return 1.
