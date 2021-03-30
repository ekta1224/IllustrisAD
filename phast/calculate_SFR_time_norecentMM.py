import numpy as np
import matplotlib
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size'] = 16
from cosmo_tools import time, snapnum2z
from plot_SFR_Z_panels_recentMM import per_bi
from matplotlib.ticker import FuncFormatter


ids = np.loadtxt('M31analogs_noMM8Gyr_mstar_noM33.txt'
#ids = np.loadtxt('M31analogs_MM1_4Gyr_mstar.txt')
print len(ids)

#remloc = ['right', 'left', 'left','left','right', 'right']#side that merger is located

for id,rem in zip(ids,remloc):
    print id
    left_sfr = []
    right_sfr = []
    times = []
    rem_sfr = []
    norem_sfr = []

    for snap in snaps:
        print int(id)
        #load gas data
        ppx2, ppy2, ppz2, ppvx2, ppvy2, ppvz2, ppm2, ppnh2, ppsfr2, ppgz2 = np.loadtxt('./SFHs/M31analog_%s_gas_properties_snap%s_rotated.txt'%(int(id), snap), usecols=(0, 1, 2, 3, 4, 5, 6, 7,8,9), unpack = True)
        print 'SFR min, max', np.min(ppsfr2), np.max(ppsfr2)

        #load orbit data
        orbit = np.loadtxt('M31analogs_major_merger_orbit_%s.txt'%int(id))
        mask = np.where(orbit[:,0] == 103)[0][0] #snap 110 = 3.96 lookback Gyr, 103= 5.08 Gyr
        xs = orbit[:,1][:mask+1]
        ys = orbit[:,2][:mask+1]
        mmt = round(time(snapnum2z(orbit[:,0][0]+1)),2)
        print mmt 

        m,b = per_bi(xs[0], ys[0], 0., 0.)
        yrange2 = m*ppx2+b
        print m,b

        print 'total lens x,y',  len(ppx2), len(ppy2)
        right2 = (ppy2 <= yrange2)
        left2 = (ppy2 > yrange2)

        print 'halves', len(ppsfr2[left2]), len(ppsfr2[right2])
        left_sfr.append(np.sum(ppsfr2[left2]))
        right_sfr.append(np.sum(ppsfr2[right2]))
        times.append(time(snapnum2z(snap)))

        if rem == 'left':
            rem_sfr.append(np.sum(ppsfr2[left2]))
            norem_sfr.append(np.sum(ppsfr2[right2]))

        if rem =='right':
            rem_sfr.append(np.sum(ppsfr2[right2]))
            norem_sfr.append(np.sum(ppsfr2[left2]))

   # np.savetxt('%s_SFR_time_summed.txt'%int(id), np.column_stack((times, left_sfr, right_sfr)), delimiter="  ")
    np.savetxt('%s_SFR_time_summed_remloc.txt'%int(id), np.column_stack((times,rem_sfr, norem_sfr)), delimiter="  ")
