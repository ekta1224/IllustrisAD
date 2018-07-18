import math
import numpy as np
from numpy.linalg import norm
from scipy.integrate import simps
from scipy.integrate import quad 

def time(z):
    '''
    Calculate lookback time for a redshift from Illustris
    '''
    t = 0.
    if z == 0.:
        t = 0.001

    if z == 1000.:
        t = 13.8

    if  0. < z < 1000.:
        h = 0.704
        H0 = h*100
        OM = 0.2726
        OL = 0.7274
        H0 = H0 * 3.241e-20 / 3.171e-17 # final units of[Gyr ^-1]

        def f(z):
            return 1/ (H0*(1+z)*np.sqrt(OM*(1+z)**3 + OL)) 

        zs = np.arange(0., z, 1e-5)
        y = f(zs)
        t = simps(y,zs)
    return t
