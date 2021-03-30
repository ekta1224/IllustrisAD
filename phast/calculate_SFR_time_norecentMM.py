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
#from plot_SFR_Z_panels_recentMM import per_bi
from matplotlib.ticker import FuncFormatter

#ids = np.loadtxt('M31analogs_noMM8Gyr_mstar.txt') #61
#print(len(ids))

ids = np.loadtxt('M31analogs_MM1_4Gyr_mstar.txt')
print(len(ids))

snaps = [134., 133., 132., 131., 129., 126., 121., 112.] # 0-4 Gyrr # , 98.] 

for id in ids:
    print(id)
    left_sfr = []
    right_sfr = []
    left_sfr_avg = []
    right_sfr_avg = []

    times = []

    for snap in snaps:
        snap = int(snap)
        #print(int(id))
        #load gas data
        ppx2, ppy2, ppz2, ppvx2, ppvy2, ppvz2, ppm2, ppnh2, ppsfr2, ppgz2 = np.loadtxt('./SFHs/M31analog_%s_gas_properties_snap%s_rotated.txt'%(int(id), snap), usecols=(0, 1, 2, 3, 4, 5, 6, 7,8,9), unpack = True)
        #print('SFR min, max', np.min(ppsfr2), np.max(ppsfr2))


        # m,b = per_bi(xs[0], ys[0], 0., 0.)
        # yrange2 = m*ppx2+b
        # print m,b

        # print 'total lens x,y',  len(ppx2), len(ppy2)
        # right2 = (ppy2 <= yrange2)
        # left2 = (ppy2 > yrange2)

        #split at x=0
        right = (ppx2 >= 0. )
        left = (ppx2 < 0.)
        
        #print('halves', len(ppsfr2[left]), len(ppsfr2[right]))
        left_sfr.append(ppsfr2[left])
        right_sfr.append(ppsfr2[right])
        times.append(time(snapnum2z(snap)))
        #print(times)
        #print(left_sfr, len(left_sfr))
        left_sfr_avg.append(np.mean(ppsfr2[left]))
        right_sfr_avg.append(np.mean(ppsfr2[right]))
    #print(left_sfr_avg, right_sfr_avg)
        # if rem == 'left':
        #     rem_sfr.append(np.sum(ppsfr2[left2]))
        #     norem_sfr.append(np.sum(ppsfr2[right2]))

        # if rem =='right':
        #     rem_sfr.append(np.sum(ppsfr2[right2]))
        #     norem_sfr.append(np.sum(ppsfr2[left2]))

    
    np.savetxt('%s_SFR_time_avg_recentMM_x0split.txt'%int(id), np.column_stack((times, left_sfr_avg, right_sfr_avg)), delimiter="  ")
