import numpy as np
from rotate_coordinates import COMdefine, RotateFrame
import os

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

    rnew_test = COMdefine(x,y,z,np.array(mass)*1e10/0.704)
    COMr.append(COMdefine(x,y,z,np.array(mass)*1e10/0.704))
    #print COMdefine(x,y,z,np.array(mass)*1e10/0.704)
    #print r1

    vnew_test = COMdefine(vx,vy,vz,np.array(mass)*1e10/0.704)
    COMv.append(COMdefine(vx,vy,vz,np.array(mass)*1e10/0.704))
    #print COMdefine(vx,vy,vz,np.array(mass)*1e10/0.704)
    #print v1

    ##### USE THIS NEW COM TO ROTATE THE PARTICLES
    m = gas[:,6]
    nh = gas[:,7]
    sfr = gas[:,8]
    gz = gas[:,9]

    x = gas[:,0] - rnew_test[0]
    y = gas[:,1] - rnew_test[1]
    z = gas[:,2] - rnew_test[2]
    r2 = np.array([x,y,z]).T 
    
    vx = gas[:,3] - vnew_test[0]
    vy = gas[:,4] - vnew_test[1]
    vz = gas[:,5] - vnew_test[2]
    v2 = np.array([vx,vy,vz]).T
        
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('../data/M31analog_%s_gas_properties_rotated_COMtest.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh, sfr, gz)), delimiter="  ")
    os.system('git add ../data/M31analog_%s_gas_properties_rotated_COMtest.txt'%id)

    #stars
    m = stars[:,6]
    tform = stars[:,7]    
    gz = stars[:,8]
    newr, newv = RotateFrame(r2, v2)

    x = stars[:,0] - rnew_test[0]
    y = stars[:,1] - rnew_test[1]
    z = stars[:,2] - rnew_test[2]
    r2 = np.array([x,y,z]).T 
    
    vx = stars[:,3] - vnew_test[0]
    vy = stars[:,4] - vnew_test[1]
    vz = stars[:,5] - vnew_test[2]
    v2 = np.array([vx,vy,vz]).T
    
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('../data/M31analog_%s_star_properties_rotated_COMtest.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, tform, gz)), delimiter="  ")
    os.system('git add ../data/M31analog_%s_star_properties_rotated_COMtest.txt'%id)


np.savetxt('COM_gas_stars_test.txt', np.column_stack((ids, COMr, COMv)))



