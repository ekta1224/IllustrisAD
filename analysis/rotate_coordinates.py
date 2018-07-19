import numpy as np

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


def COM(id):
    halo = np.loadtxt('../data/M31analogs_halo_props.txt')
    mask = np.where(halo[:,0] == id)[0][0]
    # get the halo COM position and velocity -- should this be the halo COM? or respective to each particle type?
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    
    #GAS
    gas = np.loadtxt('../data/M31analog_%s_gas_properties.txt'%id)
    x = gas[:,0] - r1[0]
    y = gas[:,1] - r1[1]
    z = gas[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = gas[:,3] - v1[0]
    vy = gas[:,4] - v1[1]
    vz = gas[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
    
    m = gas[:,6]
    nh = gas[:,7]
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('../data/M31analog_%s_gas_properties_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh)), delimiter="  ")

    #STARS
    stars = np.loadtxt('../data/M31analog_%s_star_properties.txt'%id)
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
    newr, newv = RotateFrame(r2, v2)

    np.savetxt('../data/M31analog_%s_star_properties_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, tform)), delimiter="  ")

    # DM
    dm = np.loadtxt('../data/M31analog_%s_dm_properties.txt'%id)
    x = dm[:,0] - r1[0]
    y = dm[:,1] - r1[1]
    z = dm[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = dm[:,3] - v1[0]
    vy = dm[:,4] - v1[1]
    vz = dm[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
    
    m = dm[:,6]
    newr, newv = RotateFrame(r2, v2)
    np.savetxt('../data/M31analog_%s_dm_properties_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m)), delimiter="  ")

    return 0


if __name__ == "__main__":

    COM(361428)
    #RotateFrame(r,v)
