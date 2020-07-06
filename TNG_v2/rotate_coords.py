import numpy as np
from rotate_coords_tests import RotateFrame

def rotate_original_all(id):
    '''rotates stars to angular momentum vector of stars and rotates gas
    to angular momentum vector of gas; all particles are first shifted to
    halo COM '''

    # get the halo COM position and velocity 
    halo = np.loadtxt('M31analogs_halo_props_TNG100_revised.txt')
    mask = np.where(halo[:,0] == id)[0]
    r1 = np.array([halo[:,1][mask], halo[:,2][mask], halo[:,3][mask]])
    v1 = np.array([halo[:,4][mask], halo[:,5][mask], halo[:,6][mask]])
    
    #GAS
    gas = np.loadtxt('/Users/ektapatel/Desktop/raw_particle_data/%s_gas_properties_TNGv2.txt'%id)
    m = gas[:,6]
    nh = gas[:,7]
    sfr = gas[:,8]
    gz = gas[:,9]

    #shift postions and velocities to halo COM
    x = gas[:,0] - r1[0]
    y = gas[:,1] - r1[1]
    z = gas[:,2] - r1[2]
    r2 = np.array([x,y,z]).T 
    
    vx = gas[:,3] - v1[0]
    vy = gas[:,4] - v1[1]
    vz = gas[:,5] - v1[2]
    v2 = np.array([vx,vy,vz]).T
    
    #apply rotation
    newr, newv, lnorm1 = RotateFrame(r2, v2)
    np.savetxt('/Users/ektapatel/Desktop/rotated_data/%s_gas_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, nh, sfr, gz)), delimiter="  ")


    #STARS; same procedure 
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
        
    newr, newv, lnorm2 = RotateFrame(r2, v2)
    np.savetxt('/Users/ektapatel/Desktop/rotated_data/%s_star_properties_TNGv2_rotated.txt'%id, np.column_stack((newr[:,0], newr[:,1], newr[:,2], newv[:,0], newv[:,1], newv[:,2],m, tform, gz)), delimiter="  ")
    lnormz = np.array([0,0,1])
    print(np.degrees(np.arccos(np.dot(lnorm1,lnormz))), np.degrees(np.arccos(np.dot(lnorm2,lnormz))))

    return np.degrees(np.arccos(np.dot(lnorm1,lnorm2)))

if __name__ == "__main__":

    ids = np.loadtxt('M31analogs_halo_props_TNG100_revised.txt')[:,0]
    print(ids, len(ids))
    angs = []
    for id in ids: 
        ang = rotate_original_all(int(id))
        angs.append(ang)

    np.savetxt('M31analogs_halo_angle_btwn_Lnorm_stars_gas_TNG100.txt', np.column_stack((ids, angs)))
