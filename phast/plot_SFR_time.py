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

ids = np.loadtxt('M31analogs_MM1_4Gyr_mstar.txt')


for id in ids:
    id  = int(id)
    times, left_sfr, right_sfr = np.loadtxt('%s_SFR_time_summed.txt'%id, usecols=(0,1,2), unpack=True)
    print times
    print left_sfr
    print right_sfr
    max = np.max([np.max(left_sfr), np.max(right_sfr)])

    times2 = list(times)
    times2.append(0.)
    xx = np.diff(times2)
    print xx*-1
    #widths = np.array(widths)/1e9 # these are the right #s in Gyr
    plt.figure(figsize=(8,4.5))
    ax1 = plt.subplot(211)
    ax1.bar(times, left_sfr, width=np.array(xx)*-1, align='center', alpha=0.3)#, lw=2, edgecolor='red')
    #ax1.plot(times, left_sfr, '.')
    ax1.set_xscale("log")
    ax1.set_xticks([10., 1, .1, .01])
    ax1.set_xticklabels(['10', '1', '0.1', '0.01'])
    ax1.set_ylim(0, max+1)
    ax1.set_xlim(ax1.get_xlim()[::-1]) 
    ax1.set_title('left')
    #ax1.set_xlabel('Age [Gyr]', labelpad=-3)
    ax1.set_ylabel(r'SFR [M$_{\odot}$/yr]') 
    #plt.savefig('recentMM_SFR_time_left.pdf')

    #plt.figure(figsize=(8,4.5))
    ax2 = plt.subplot(212)
    ax2.bar(times, right_sfr, width=np.array(xx)*-1, align='center', color='red', alpha=0.3)
    #ax2.plot(times, right_sfr, '.')
    ax2.set_xscale("log")
    ax2.set_xticks([10., 1, .1, .01])
    ax2.set_xticklabels(['10', '1', '0.1', '0.01'])
    ax2.set_xlim(ax2.get_xlim()[::-1]) 
    ax2.set_ylim(0, max+1)
    ax2.set_title('right')
    ax2.set_xlabel('Age [Gyr]', labelpad=-3)
    ax2.set_ylabel(r'SFR [M$_{\odot}$/yr]') 
    plt.tight_layout()
    plt.savefig('recentMM_SFR_time_%s.pdf'%int(id))
