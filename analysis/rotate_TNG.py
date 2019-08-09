import numpy as np 

#Ekta's functions to calculate the CoM and to rotate the frame -- as in /Volumes/FRIEND/analogs/IllustrisAD/analysis/rotate_coordinates.py
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

#read in IDs
isolated_ids = np.loadtxt('/Users/amandaquirk/Desktop/isolated_central_IDs.txt')
M33_comp_ids = np.loadtxt('/Users/amandaquirk/Desktop/M33com_central_IDs.txt')
all_ids = list(isolated_ids) + list(M33_comp_ids)

star_x, star_y, star_z, star_vx, star_vy, star_vz, star_tform = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_star_properties.txt'.format(all_ids[0]), usecols=(0,1,2,3,4,5,7), unpack=True)
gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_nf = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_gas_properties.txt'.format(all_ids[0]), usecols=(0,1,2,3,4,5,7), unpack=True)
dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz, dm_potential = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_dm_properties.txt'.format(all_ids[0]), usecols=(0,1,2,3,4,5), unpack=True)

def potential_to_mass(potential):
	return 

dm_mass = potential_to_mass(dm_potential)

#step 1 is to find the CoM of the dark matter -- set this as the origin
vx_com, vy_com, vz_com = COMdefine(dm_vx, dm_vy, dm_vz, dm_mass)
x_com, y_com, z_com = COMdefine(dm_x, dm_y, dm_z, dm_mass)
#print('Found COMs')

#subtract off COM so all particles are relative to that origin
star_x0 = star_x - x_com
star_y0 = star_y - y_com
star_z0 = star_z - z_com
star_vx0 = star_vx - vx_com
star_vy0 = star_vy - vy_com
star_vz0 = star_vz - vz_com
gas_x0 = gas_x - x_com
gas_y0 = gas_y - y_com
gas_z0 = gas_z - z_com
gas_vx0 = gas_vx - vx_com
gas_vy0 = gas_vy - vy_com
gas_vz0 = gas_vz - vz_com
#print('Shifted coordinates to star COM')

r_star = np.array([star_x0, star_y0, star_z0]).T #array to go into RotateFrame
v_star = np.array([star_vx0, star_vy0, star_vz0]).T #array to go into RotateFrame
r_gas = np.array([gas_x0, gas_y0, gas_z0]).T #array to go into RotateFrame
v_gas = np.array([gas_vx0, gas_vy0, gas_vz0]).T #array to go into RotateFrame

#transform the coordinates
star_r_transf, star_v_transf = RotateFrame(r_star, v_star)
#print('Rotated star coordinates')
gas_r_transf, gas_v_transf = RotateFrame(r_gas, v_gas)
#print('Rotated gas coordinates')

# #save gas and star data
np.savetxt('/Volumes/FRIEND/analogs/TNGdata/{}_rotated_star_particles.txt'.format(int(all_ids[0])), np.column_stack([star_r_transf[:,0], star_r_transf[:,1], star_r_transf[:,2], star_v_transf[:,0], star_v_transf[:,1], star_v_transf[:,2], star_tform]))
np.savetxt('/Volumes/FRIEND/analogs/TNGdata/{}_rotated_gas_particles.txt'.format(int(all_ids[0])), np.column_stack([gas_r_transf[:,0], gas_r_transf[:,1], gas_r_transf[:,2], gas_v_transf[:,0], gas_v_transf[:,1], gas_v_transf[:,2], gas_nf]))




