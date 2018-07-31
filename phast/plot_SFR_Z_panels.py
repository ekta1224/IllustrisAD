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
#matplotlib.rcParams['text.usetex'] = True

#IDs for M31 analogs with not major merger in the last 8 Gyr and a stellar mass cut of 5e10-2e11 applied
ids = np.loadtxt('M31analogs_noMM8Gyr_mstar.txt')

left_sfr = []
right_sfr = []
left_z = []
right_z = []


for id in ids:
    id = int(id)
    print id
    #load stellar data
    ppx, ppy, ppz, ppvx, ppvy, ppvz, ppm, ppage, ppgz = np.loadtxt('../data/M31analog_%s_star_properties_rotated.txt'%id, usecols=(0, 1, 2, 3, 4, 5, 6, 7,8), unpack = True)
    #print len(ppz)
    #print len(ppgz), len(ppage[ppage >0.])

    ppx = ppx[ppage > 0.]
    ppy = ppy[ppage > 0.]
    ppz = ppz[ppage > 0.]
    ppm = ppm[ppage > 0.]
    ppgz = np.abs(ppgz[ppage > 0.])
    ppage = ppage[ppage > 0.]

    #load gas data
    ppx2, ppy2, ppz2, ppvx2, ppvy2, ppvz2, ppm2, ppnh2, ppsfr2, ppgz2 = np.loadtxt('../data/M31analog_%s_gas_properties_rotated.txt'%id, usecols=(0, 1, 2, 3, 4, 5, 6, 7,8,9), unpack = True)
    print np.min(ppsfr2), np.max(ppsfr2)

    # plot
    plt.figure(figsize=(13,6))
    plt.subplot(121)
    plt.scatter(ppx, ppy, c=ppgz, marker='.', s=1,vmin=0, vmax=np.max(ppgz), cmap=plt.cm.get_cmap('viridis_r', 8))#,alpha=0.5)
    plt.plot(0., 0., 'kx',ms=8)
    plt.axvline(x=0, ls='--', color='gray')
    plt.colorbar(label=r'metallicity $M_Z/M_{gas}$')
    plt.xlim(-25,25)
    plt.ylim(-25,25)
    plt.title('%s (stars)'%id)
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')


    plt.subplot(122)
    plt.scatter(ppx2, ppy2, c=ppsfr2, marker='.', s=2,norm=colors.LogNorm(vmin=1e-4, vmax=np.max(ppsfr2)), cmap='viridis_r')
    plt.plot(0., 0., 'kx',ms=8)
    plt.axvline(x=0, ls='--', color='gray')
    plt.colorbar(label=r'SFR [M$_{\odot}$/yr]')
    plt.xlim(-25,25)
    plt.ylim(-25,25)
    plt.title('%s (gas)'%id)
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')
    plt.tight_layout()
    #plt.savefig('%s_sfr_z.pdf'%id)
    plt.savefig('./norecentMMplots/%s_sfr_z.png'%id, dpi=300)
    plt.close()

    #save the mean SFR and metallicity of each half
    print len(ppx), len(ppy)
    right = (ppx >= 0. )
    left = (ppx < 0.)

    print len(ppx2), len(ppy2)
    right2 = (ppx2 >= 0.)
    left2 = (ppx2 < 0.)
    
    print len(ppgz[left]), len(ppgz[right])
    print len(ppsfr2[left2]), len(ppsfr2[right2])

    left_z.append(np.mean(ppgz[left]))
    right_z.append(np.mean(ppgz[right]))
    left_sfr.append(np.mean(ppsfr2[left2]))
    right_sfr.append(np.mean(ppsfr2[right2]))

np.savetxt('M31analogs_norecentMM_mean_props.txt', np.column_stack((left_z, right_z, left_sfr, right_sfr)), delimiter="  ")
    
    
