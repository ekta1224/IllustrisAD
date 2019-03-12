import numpy as np
from rotate_coordinates import COMdefine
import matplotlib.pyplot as plt
from adjustText import adjust_text

#check the diffference between COM calculations using DM, star, and
#gas particles for a subset of halos

ids = [356427, 365141, 383739, 385230, 385343, 387913, 388958, 389486, 390669, 392084, 394566, 395868, 396624, 397062, 397148, 397204, 397889, 398622, 400733, 401663]

halo = np.loadtxt('../data/M31analogs_halo_props.txt')
dx = []
dy = []
dz = []

dvx = []
dvy = []
dvz = []

dx2 = []
dy2 = []
dz2 = []

dvx2 = []
dvy2 = []
dvz2 = []

dx3 = []
dy3 = []
dz3 = []

dvx3 = []
dvy3 = []
dvz3 = []


for id in ids:
    print id
    mask = np.where(halo[:,0] == id)[0][0]

    #halo COM
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    #print r1

    #load gas
    gas = np.loadtxt('../data/M31analog_%s_gas_properties.txt'%id)
    m = gas[:,6]

    gasr1 = COMdefine(gas[:,0], gas[:,1], gas[:,2], m*1e10/0.704)
    gasv1 = COMdefine(gas[:,3], gas[:,4], gas[:,5], m*1e10/0.704)
    #print gasr1

    #load stars
    stars = np.loadtxt('../data/M31analog_%s_star_properties.txt'%id)
    m = stars[:,6]

    starsr1 = COMdefine(stars[:,0], stars[:,1], stars[:,2], m*1e10/0.704)
    starsv1 = COMdefine(stars[:,3], stars[:,4], stars[:,5], m*1e10/0.704)
    #print starsr1


    #print v1
    #print gasv1
    #print starsv1

    dx.append(np.abs(np.abs(gasr1[0])-np.abs(r1[0])))
    dy.append(np.abs(np.abs(gasr1[1])-np.abs(r1[1])))
    dz.append(np.abs(np.abs(gasr1[2])-np.abs(r1[2])))

    dvx.append(np.abs(np.abs(gasv1[0])-np.abs(v1[0])))
    dvy.append(np.abs(np.abs(gasv1[1])-np.abs(v1[1])))
    dvz.append(np.abs(np.abs(gasv1[2])-np.abs(v1[2])))

    dx2.append(np.abs(np.abs(starsr1[0])-np.abs(r1[0])))
    dy2.append(np.abs(np.abs(starsr1[1])-np.abs(r1[1])))
    dz2.append(np.abs(np.abs(starsr1[2])-np.abs(r1[2])))

    dvx2.append(np.abs(np.abs(starsv1[0])-np.abs(v1[0])))
    dvy2.append(np.abs(np.abs(starsv1[1])-np.abs(v1[1])))
    dvz2.append(np.abs(np.abs(starsv1[2])-np.abs(v1[2])))

    dx3.append(np.abs(np.abs(starsr1[0])-np.abs(gasr1[0])))
    dy3.append(np.abs(np.abs(starsr1[1])-np.abs(gasr1[1])))
    dz3.append(np.abs(np.abs(starsr1[2])-np.abs(gasr1[2])))

    dvx3.append(np.abs(np.abs(starsv1[0])-np.abs(gasv1[0])))
    dvy3.append(np.abs(np.abs(starsv1[1])-np.abs(gasv1[1])))
    dvz3.append(np.abs(np.abs(starsv1[2])-np.abs(gasv1[2])))



### SCATTER PLOTS ###
plt.figure()
plt.subplot(121)
plt.scatter(dx, dy, c=dz, cmap='magma')
plt.xlabel(r'$\Delta$ X [kpc]')
plt.ylabel(r'$\Delta$ Y [kpc]')
plt.colorbar(label=r'$\Delta$ Z [kpc]')
texts = [plt.text(dx[i], dy[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)

plt.subplot(122)
plt.scatter(dvx, dvy, c=dvz, cmap='viridis')
plt.xlabel(r'$\Delta V_x$ [km/s]')
plt.ylabel(r'$\Delta V_y$ [km/s]')
plt.colorbar(label=r'$\Delta V_z$ [km/s]')
texts = [plt.text(dvx[i], dvy[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)
plt.suptitle('COM difference using gas particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_gas_scatter.pdf')
plt.close()

## stars
plt.figure()
plt.subplot(121)
plt.scatter(dx2, dy2, c=dz2, cmap='magma')
plt.xlabel(r'$\Delta$ X [kpc]')
plt.ylabel(r'$\Delta$ Y [kpc]')
plt.colorbar(label=r'$\Delta$ Z [kpc]')
texts = [plt.text(dx2[i], dy2[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)

plt.subplot(122)
plt.scatter(dvx2, dvy2, c=dvz2, cmap='viridis')
plt.xlabel(r'$\Delta V_x$ [km/s]')
plt.ylabel(r'$\Delta V_y$ [km/s]')
plt.colorbar(label=r'$\Delta V_z$ [km/s]')
texts = [plt.text(dvx2[i], dvy2[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)
plt.suptitle('COM difference using star particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_stars_scatter.pdf')
plt.close()

##stars-gas
plt.figure()
plt.subplot(121)
plt.scatter(dx3, dy3, c=dz3, cmap='magma')
plt.xlabel(r'$\Delta$ X [kpc]')
plt.ylabel(r'$\Delta$ Y [kpc]')
plt.colorbar(label=r'$\Delta$ Z [kpc]')
texts = [plt.text(dx3[i], dy3[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)

plt.subplot(122)
plt.scatter(dvx3, dvy3, c=dvz3, cmap='viridis')
plt.xlabel(r'$\Delta V_x$ [km/s]')
plt.ylabel(r'$\Delta V_y$ [km/s]')
plt.colorbar(label=r'$\Delta V_z$ [km/s]')
texts = [plt.text(dvx3[i], dvy3[i], '%s'%ids[i], fontsize=7.5) for i in range(len(ids))]
adjust_text(texts)
plt.suptitle('COM difference btwn gas and star particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_stars_gas_scatter.pdf')
plt.close()



### HISTOGRAMS ###

plt.figure()
plt.subplot(231)
plt.hist(dx, bins=5)
plt.xlabel(r'$\Delta$ X [kpc]')

plt.subplot(232)
plt.hist(dy, bins=5)
plt.xlabel(r'$\Delta$ Y [kpc]')

plt.subplot(233)
plt.hist(dz, bins=5)
plt.xlabel(r'$\Delta$ Z [kpc]')

plt.subplot(234)
plt.hist(dvx, bins=5)
plt.xlabel(r'$\Delta V_x$  [km/s]')

plt.subplot(235)
plt.hist(dvy, bins=5)
plt.xlabel(r'$\Delta V_y$  [km/s]')

plt.subplot(236)
plt.hist(dvz, bins=5)
plt.xlabel(r'$\Delta V_z$  [km/s]')

plt.suptitle('COM difference using gas particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_gas.pdf')
plt.close()


## repeat with stars
plt.figure()
plt.subplot(231)
plt.hist(dx2, bins=5)
plt.xlabel(r'$\Delta$ X [kpc]')

plt.subplot(232)
plt.hist(dy2, bins=5)
plt.xlabel(r'$\Delta$ Y [kpc]')

plt.subplot(233)
plt.hist(dz2, bins=5)
plt.xlabel(r'$\Delta$ Z [kpc]')

plt.subplot(234)
plt.hist(dvx2, bins=5)
plt.xlabel(r'$\Delta V_x$  [km/s]')

plt.subplot(235)
plt.hist(dvy2, bins=5)
plt.xlabel(r'$\Delta V_y$  [km/s]')

plt.subplot(236)
plt.hist(dvz2, bins=5)
plt.xlabel(r'$\Delta V_z$  [km/s]')

plt.suptitle('COM difference using star particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_stars.pdf')

## stars-gas
plt.figure()
plt.subplot(231)
plt.hist(dx3, bins=5)
plt.xlabel(r'$\Delta$ X [kpc]')

plt.subplot(232)
plt.hist(dy3, bins=5)
plt.xlabel(r'$\Delta$ Y [kpc]')

plt.subplot(233)
plt.hist(dz3, bins=5)
plt.xlabel(r'$\Delta$ Z [kpc]')

plt.subplot(234)
plt.hist(dvx3, bins=5)
plt.xlabel(r'$\Delta V_x$  [km/s]')

plt.subplot(235)
plt.hist(dvy3, bins=5)
plt.xlabel(r'$\Delta V_y$  [km/s]')

plt.subplot(236)
plt.hist(dvz3, bins=5)
plt.xlabel(r'$\Delta V_z$  [km/s]')

plt.suptitle('COM difference btwn gas and star particles')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('COM_test_stars_gas.pdf')


