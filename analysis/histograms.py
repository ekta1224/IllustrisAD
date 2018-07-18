import numpy as np
import matplotlib.pyplot as plt
from redshift2time import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def nh_hist(filename):
    ''' plot histograms for the netural hydrogen fraction in gas cells

    input: full filename
    '''
    gas = np.loadtxt('../data/%s'%filename)
    nh = gas[:,7]
    print 'number of gas cells', len(nh)
    plt.figure()
    plt.hist(nh, bins=20, histtype='step', lw=2, weights=np.ones_like(nh)/len(nh))
    plt.xlabel('neutral hydrogren fraction')
    plt.title('M31 analog - Illustris %s'%filename[10:15])
    plt.savefig('nh_hist_%s.pdf'%filename[10:15])
    plt.close()
    return 0

def age_hist(filename):
    
    ''' plot histograms for stellar age (taken from stellar formation time)
    '''
    stars = np.loadtxt('../data/%s'%filename)
    aform = np.array(stars[:,7]) #a= 1/(1+z) --> z = 1/a - 1
    aform = aform[aform >= 0.]
    print 'number of star particles', len(aform)
    zform = [(1./a) -1. for a in aform]
    cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
    tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
    tage = np.array([13.8-t for t in tform])
    plt.figure()
    plt.hist(tage, bins=20, histtype='step', lw=2, weights=np.ones_like(tage)/len(tage), color='crimson')
    plt.xlabel('stellare ages [Gyr]')
    plt.title('M31 analog - Illustris %s'%filename[10:15])
    plt.savefig('stellar_age_hist_%s.pdf'%filename[10:15])
    plt.close()

    return 0
             
if __name__ == "__main__":
    filename = 'M31analog_361428_star_properties.txt'
    age_hist(filename)

    filename = 'M31analog_361428_gas_properties.txt'
    nh_hist(filename)

    #break sets into ages of 30 Myr, 400 Myr, 2 Gyr, 4 Gyr?
