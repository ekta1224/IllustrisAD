import numpy as np

data = np.loadtxt('../data/M31analogs_halo_props.txt')
ids = data[:,0]
mmt = data[:,9]
mstar = data[:,10]/0.704*1e10

print np.min(mstar), np.max(mstar)

print len(mmt)
print len(mmt[(mmt >=1.)*(mmt <=4.)*(mstar >= 5e10)*(mstar <= 2e11)])
print len(mmt[(mmt >= 8.)*(mstar >= 5e10)*(mstar <= 2e11)])

np.savetxt('M31analogs_noMM8Gyr_mstar.txt', ids[(mmt >= 8.)*(mstar >= 5e10)*(mstar <= 2e11)])
