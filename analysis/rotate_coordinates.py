import numpy as np
import os

def COMdefine(a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)
    
        # xcomponent Center of mass  
        Acom = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.sum(b*m)/np.sum(m)
        # zcomponent
        Ccom = np.sum(c*m)/np.sum(m)

        return Acom, Bcom, Ccom

def RotateFrame(posI,velI):
    # input:  3D array of positions and velocities
    # returns: 3D array of rotated positions and velocities such that j is in z direction

    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel


def COM(id, noM33=False):
    halo = np.loadtxt('../data/COM_gas_stars_test.txt')#M31analogs_halo_props.txt')
    if noM33:
	    halo = np.loadtxt('../data/M31analogs_halo_props_strictly_noM33.txt')
    mask = np.where(halo[:,0] == id)[0]#[0]
    # get the halo COM position and velocity -- should this be the halo COM? or respective to each particle type?
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    #print r1, v1

    #GAS
    gas = np.loadtxt('../data/M31analog_%s_gas_properties.txt'%id)
    m = gas[:,6]
    nh = gas[:,7]
    sfr = gas[:,8]
    gz = gas[:,9]

    x = gas[:,0] - r1[0]
    y = gas[:,1] - r1[1]
    z = gas[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = gas[:,3] - v1[0]
    vy = gas[:,4] - v1[1]
    vz = gas[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
        
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('/Volumes/Titan/analogs/data/COM_test/M31analog_%s_gas_properties_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh, sfr, gz)), delimiter="  ")
    

    #STARS
    stars = np.loadtxt('/Volumes/Titan/analogs/data/COM_test/M31analog_%s_star_properties.txt'%id)
    x = stars[:,0] - r1[0]
    y = stars[:,1] - r1[1]
    z = stars[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = stars[:,3] - v1[0]
    vy = stars[:,4] - v1[1]
    vz = stars[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
    
    m = stars[:,6]
    tform = stars[:,7]    
    gz = stars[:,8]
    newr, newv = RotateFrame(r2, v2)

    np.savetxt('../data/M31analog_%s_star_properties_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, tform, gz)), delimiter="  ")

    return 0

# def COM_calculate(id, snap):
#     ''' this is to rotate the data obtained to create SFHs only'''

#     #GAS
#     gas = np.loadtxt('../phast/SFHs/M31analog_%s_gas_properties_snap%s.txt'%(id,snap))
#     m = gas[:,6]
#     nh = gas[:,7]
#     sfr = gas[:,8]
#     gz = gas[:,9]

#     r1 = COMdefine(gas[:,0], gas[:,1], gas[:,2], m*1e10/0.704)
#     v1 = COMdefine(gas[:,3], gas[:,4], gas[:,5], m*1e10/0.704)

#     x = gas[:,0] - r1[0]
#     y = gas[:,1] - r1[1]
#     z = gas[:,2] - r1[2]
#     r2 = np.array([x,y,z]).T 
    
#     vx = gas[:,3] - v1[0]
#     vy = gas[:,4] - v1[1]
#     vz = gas[:,5] - v1[2]
#     v2 = np.array([vx,vy,vz]).T
        
#     newr, newv = RotateFrame(r2, v2)
#     np.savetxt('../phast/SFHs/M31analog_%s_gas_properties_snap%s_rotated.txt'%(id, snap), np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh, sfr, gz)), delimiter="  ")
    
#     return 0

# if __name__ == "__main__":
#     ids = np.loadtxt('../phast/M31analogs_MM1_4Gyr_mstar.txt')
#     print len(ids)
#     snaps = [98, 112, 121, 126, 129, 131, 132, 133, 134]
#     for id in ids:
#         for snap in snaps:
#             COM_calculate(int(id), snap)
#     assert False

#     #rotate additional analogs that don't have to have an M33
#     orig_subids = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt')
#     add_subids = np.concatenate((np.loadtxt('../phast/M31analogs_MM1_4Gyr_mstar_noM33.txt'),np.loadtxt('../phast/M31analogs_noMM8Gyr_mstar_noM33.txt')))
#     print len(add_subids)
#     subids = []
#     for i in add_subids:
#         if i in orig_subids:
#             continue
#         else:
#             subids.append(i)
#     print len(subids)

#     for id in subids:
#         print id
#         COM(int(id), noM33=True)

#     assert False
#     ids = np.loadtxt('../data/M31analogs_halo_props.txt')[:,0]
#     for id in ids:
#         print id
#         COM(int(id))

#run above functions
halos = [3.564270000000000000e+05, 3.874560000000000000e+05, 4.000040000000000000e+05, 4.053000000000000000e+05, 4.089160000000000000e+05, 4.122090000000000000e+05, 4.174200000000000000e+05, 4.236530000000000000e+05, 4.282060000000000000e+05, 4.341020000000000000e+05]
halos = [int(a) for a in halos]

for halo in halos: 
    COM(halo, noM33=True) 