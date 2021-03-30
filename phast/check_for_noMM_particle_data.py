import numpy as np
import os

ids = np.loadtxt('M31analogs_noMM8Gyr_mstar_noM33.txt')
for id in ids:
    id = int(id)
    print(id)
    os.system('cp ../data/M31analog_%s_gas_properties_rotated.txt ./SFHs'%id)
#    ppx2, ppy2, ppz2, ppvx2, ppvy2, ppvz2, ppm2, ppnh2, ppsfr2, ppgz2 = np.loadtxt('../data/M31analog_%s_gas_properties_rotated.txt'%id,usecols=(0, 1, 2, 3, 4, 5, 6, 7,8,9), unpack = True)
