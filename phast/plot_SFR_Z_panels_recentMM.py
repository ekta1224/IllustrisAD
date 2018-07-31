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

def per_bi(x1, y1, x2, y2):
    '''return the slope and intercept of the perpendicular bisector between two points'''
    mid = ((x1+x2)/2., (y1+y2)/2.)
    print mid
    slope = (y2-y1)/(x2-x1)
    print slope
    m = (1./slope)*-1.
    b = mid[1] - (mid[0]*m)
    return m,b

    
#IDs for M31 analogs with a major merger in the last 1-4 Gyr and a stellar mass cut of 5e10-2e11 applied
ids = np.loadtxt('M31analogs_MM1_4Gyr_mstar.txt')

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

    #load orbit data
    orbit = np.loadtxt('M31analogs_major_merger_orbit_%s.txt'%id)
    mask = np.where(orbit[:,0] == 103)[0][0] #snap 110 = 3.96 lookback Gyr, 103= 5.08 Gyr
    xs = orbit[:,1][:mask+1]
    ys = orbit[:,2][:mask+1]

    # plot
    plt.figure(figsize=(13,6))
    plt.subplot(121)
    plt.scatter(ppx, ppy, c=ppgz, marker='.', s=1,vmin=0, vmax=np.max(ppgz), cmap=plt.cm.get_cmap('viridis_r', 8))#,alpha=0.5)
    plt.plot(0., 0., 'kx',ms=8)
    plt.plot(xs ,ys ,'k--')
    plt.plot(xs[0], ys[0], 'ko')
    m,b = per_bi(xs[0], ys[0], 0., 0.)
    yrange = m*ppx+b
    plt.plot(ppx, yrange, '--', c='gray')
    plt.colorbar(label=r'metallicity $M_Z/M_{gas}$')
    plt.xlim(-75,75)
    plt.ylim(-75,75)
    plt.title('%s (stars)'%id)
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')
                   
    plt.subplot(122)
    plt.scatter(ppx2, ppy2, c=ppsfr2, marker='.', s=2,norm=colors.LogNorm(vmin=1e-4, vmax=np.max(ppsfr2)), cmap='viridis_r')
    plt.plot(0., 0., 'kx', ms=8)
    plt.plot(xs , ys ,'k--')
    plt.plot(xs[0], ys[0], 'ko')
    yrange2 = m*ppx2+b
    plt.plot(ppx2, yrange2, '--', c='gray')
    plt.colorbar(label=r'SFR [M$_{\odot}$/yr]')
    plt.xlim(-75,75)
    plt.ylim(-75,75)
    plt.title('%s (gas)'%id)
    plt.xlabel('X [kpc]')
    plt.ylabel('Y [kpc]')
    plt.tight_layout()
    plt.figtext(0.6, 0.85, 'merger time: %s Gyr' %round(time(snapnum2z(orbit[:,0][0]+1)),2))
    #plt.savefig('%s_sfr_z.pdf'%id)
    #plt.savefig('./norecentMMplots/%s_sfr_z.png'%id, dpi=300)
    plt.savefig('./recentMMplots/%s_recentMM_sfr_z.png'%id, dpi=300)
    plt.close()

    #save the mean SFR and metallicity of each half
    print len(ppx), len(ppy)
    right = (ppy <= yrange)
    left = (ppy > yrange)

    print len(ppx2), len(ppy2)
    right2 = (ppy2 <= yrange2)
    left2 = (ppy2 > yrange2)
    
    print len(ppgz[left]), len(ppgz[right])
    print len(ppsfr2[left2]), len(ppsfr2[right2])

    left_z.append(np.mean(ppgz[left]))
    right_z.append(np.mean(ppgz[right]))
    left_sfr.append(np.mean(ppsfr2[left2]))
    right_sfr.append(np.mean(ppsfr2[right2]))

np.savetxt('M31analogs_recentMM_mean_props.txt', np.column_stack((left_z, right_z, left_sfr, right_sfr)), delimiter="  ")
