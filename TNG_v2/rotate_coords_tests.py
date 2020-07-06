import numpy as np

#original function to rotate frame, applied to stellar and gaseous data individually 

def RotateFrame(posI,velI):
    # input:  3D array of positions and velocities
    # returns: 3D array of rotated positions and velocities such that j is in z direction

    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))
    #print(L_norm)

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
    
    return pos, vel, L_norm


def rotate_original(id):
    halo = np.loadtxt('M31analogs_halo_props_TNG100_revised.txt')
    mask = np.where(halo[:,0] == id)[0]
    # get the halo COM position and velocity 
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    
    #GAS
    gas = np.loadtxt('/Users/ektapatel/Desktop/raw_particle_data/%s_gas_properties_TNGv2.txt'%id)
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
    np.savetxt('/Users/ektapatel/Desktop/coord_rotation_tests/original_rotation/%s_gas_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh, sfr, gz)), delimiter="  ")


    #STARS
    stars = np.loadtxt('/Users/ektapatel/Desktop/raw_particle_data/%s_star_properties_TNGv2.txt'%id)
    m = stars[:,6]
    tform = stars[:,7]    
    gz = stars[:,8]

    x = stars[:,0] - r1[0]
    y = stars[:,1] - r1[1]
    z = stars[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = stars[:,3] - v1[0]
    vy = stars[:,4] - v1[1]
    vz = stars[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
        
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('/Users/ektapatel/Desktop/coord_rotation_tests/original_rotation/%s_star_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, tform, gz)), delimiter="  ")

    return 0

def rotate_modified(id):
    #rotate particles/cells using the combined angular momentum of the stars and gas
    halo = np.loadtxt('M31analogs_halo_props_TNG100_revised.txt')
    mask = np.where(halo[:,0] == id)[0]
    # get the halo COM position and velocity 
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    
    #GAS
    gas = np.loadtxt('/Users/ektapatel/Desktop/raw_particle_data/%s_gas_properties_TNGv2.txt'%id)
    m_gas = gas[:,6]
    nh_gas = gas[:,7]
    sfr_gas = gas[:,8]
    gz_gas = gas[:,9]

    xg = gas[:,0] - r1[0]
    yg = gas[:,1] - r1[1]
    zg = gas[:,2] - r1[2]
    
    vxg = gas[:,3] - v1[0]
    vyg = gas[:,4] - v1[1]
    vzg = gas[:,5] - v1[2]

        
    #newr, newv = RotateFrame(r2, v2)

    #STARS
    stars = np.loadtxt('/Users/ektapatel/Desktop/raw_particle_data/%s_star_properties_TNGv2.txt'%id)
    m = stars[:,6]
    tform = stars[:,7]    
    gz = stars[:,8]

    x = stars[:,0] - r1[0]
    y = stars[:,1] - r1[1]
    z = stars[:,2] - r1[2]
    
    vx = stars[:,3] - v1[0]
    vy = stars[:,4] - v1[1]
    vz = stars[:,5] - v1[2]
    
        
    gas_ind = len(xg)
    star_ind = len(x)
    x_tot = list(xg) + list(x)
    y_tot = list(yg) + list(y)
    z_tot = list(zg) + list(z)

    vx_tot = list(vxg) + list(vx)
    vy_tot = list(vyg) + list(vy)
    vz_tot = list(vzg) + list(vz)

    r2 = np.array([x_tot,y_tot,z_tot]).T 
    v2 = np.array([vx_tot,vy_tot,vz_tot]).T

    newr, newv = RotateFrame(r2, v2)
    print(len(newr[:,0][:gas_ind]), len(gz_gas))
    print(len(newr[:,0][:star_ind]), len(gz))
    np.savetxt('/Users/ektapatel/Desktop/coord_rotation_tests/modified_rotation/%s_gas_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0][:gas_ind], newr[:,1][:gas_ind], newr[:,2][:gas_ind], newv[:,0][:gas_ind], newv[:,1][:gas_ind], newv[:,2][:gas_ind],m_gas, nh_gas, sfr_gas, gz_gas)), delimiter="  ")

    np.savetxt('/Users/ektapatel/Desktop/coord_rotation_tests/modified_rotation/%s_star_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0][:star_ind], newr[:,1][:star_ind], newr[:,2][:star_ind], newv[:,0][:star_ind], newv[:,1][:star_ind], newv[:,2][:star_ind],m, tform, gz)), delimiter="  ")


    return 0
if __name__ == "__main__":

    test_ids = [455957, 466746, 478298, 487042, 499113, 479264, 467575, 458864, 489593, 504142] 
    for id in test_ids:
        rotate_original(id)
        rotate_modified(id)
        
    
