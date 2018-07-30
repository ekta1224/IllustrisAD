import numpy as np

data = np.loadtxt('../data/M31analogs_halo_props.txt')
mmt = data[:,9]
print mmt

print len(mmt)
print len(mmt[(mmt >=1.)*(mmt <=4.)])
print len(mmt[(mmt >= 8.)]), mmt[mmt >= 8.]
