import numpy as np
from rotate_coordinates import COMdefine


#calculte the combined gas and stellar COM for some subset of analogs
#and redo the full analysis to see whether the rotation curves change

#subset of 10 from the primary sample
ids = [3.564270000000000000e+05,3.874560000000000000e+05,4.000040000000000000e+05,4.053000000000000000e+05,4.089160000000000000e+05,4.122090000000000000e+05,4.174200000000000000e+05,4.236530000000000000e+05,4.282060000000000000e+05,4.341020000000000000e+05]

halo = np.loadtxt('../data/M31analogs_halo_props.txt')

COMr = []
COMv = []
for id in ids:
    print(id)
    id = int(id)
    mask = np.where(halo[:,0] == id)[0][0]

    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])

    #load gas
    gas = np.loadtxt('../data/M31analog_%s_gas_properties.txt'%id)
    mgas = gas[:,6]

    #gasr1 = COMdefine(gas[:,0], gas[:,1], gas[:,2], m*1e10/0.704)
    #gasv1 = COMdefine(gas[:,3], gas[:,4], gas[:,5], m*1e10/0.704)
    

    #load stars
    stars = np.loadtxt('../data/M31analog_%s_star_properties.txt'%id)
    mstars = stars[:,6]

    #starsr1 = COMdefine(stars[:,0], stars[:,1], stars[:,2], m*1e10/0.704)
    #starsv1 = COMdefine(stars[:,3], stars[:,4], stars[:,5], m*1e10/0.704)

    print(len(gas[:,0]), len(stars[:,0]))
    x = list(gas[:,0]) + list(stars[:,0])
    y = list(gas[:,1]) + list(stars[:,1])
    z = list(gas[:,2]) + list(stars[:,2])

    vx = list(gas[:,3]) + list(stars[:,3])
    vy = list(gas[:,4]) + list(stars[:,4])
    vz = list(gas[:,5]) + list(stars[:,5])

    mass = list(mgas) + list(mstars)

    COMr.append(COMdefine(x,y,z,np.array(mass)*1e10/0.704))
    #print COMdefine(x,y,z,np.array(mass)*1e10/0.704)
    #print r1

    COMv.append(COMdefine(vx,vy,vz,np.array(mass)*1e10/0.704))
    #print COMdefine(vx,vy,vz,np.array(mass)*1e10/0.704)
    #print v1

np.savetxt('COM_gas_stars_test.txt', np.column_stack((ids, COMr, COMv)))
