import numpy as np 
import matplotlib.pyplot as plt 

#Ekta's function to rotate the frame -- as in /Volumes/FRIEND/analogs/IllustrisAD/analysis/rotate_coordinates.py
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
noM33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_noM33_TNG100.txt')
M33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_M33_TNG100.txt')
all_ids = list(noM33_ids)# + list(M33_ids)

#read in all of the subhalo properties to get the COMs
ref_ID, x_coms, y_coms, z_coms, vx_coms, vy_coms, vz_coms = np.loadtxt('/Users/amandaquirk/Desktop/M31analogs_halo_props_TNG100.txt', usecols=(0,1,2,3,4,5,6,), unpack=True)

for i in range(len(all_ids)):
	star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass, star_tform = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_star_properties.txt'.format(all_ids[i]), usecols=(0,1,2,3,4,5,6,7), unpack=True)
	gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_nf = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_gas_properties.txt'.format(all_ids[i]), usecols=(0,1,2,3,4,5,7), unpack=True)
	#dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/{}_dm_properties.txt'.format(all_ids[i]), usecols=(0,1,2,3,4,5), unpack=True)

	#dm_mass = np.ones(len(dm_x)) * 7.5 * 10**6 #constant throughout the simulation

	#grab the subhalo's position and velocity COM (really the position but we're going with it)
	N = np.where(all_ids[i] == ref_ID)
	#print(N)
	x_com = x_coms[N]
	y_com = y_coms[N]
	z_com = z_coms[N]
	vx_com = vx_coms[N]
	vy_com = vy_coms[N]
	vz_com = vz_coms[N]

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
	#print('Shifted coordinates to COM')

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
	np.savetxt('/Volumes/FRIEND/analogs/TNGdata/{}_rotated_star_COM_particles.txt'.format(int(all_ids[i])), np.column_stack([star_r_transf[:,0], star_r_transf[:,1], star_r_transf[:,2], star_v_transf[:,0], star_v_transf[:,1], star_v_transf[:,2], star_tform]))
	np.savetxt('/Volumes/FRIEND/analogs/TNGdata/{}_rotated_gas_COM_particles.txt'.format(int(all_ids[i])), np.column_stack([gas_r_transf[:,0], gas_r_transf[:,1], gas_r_transf[:,2], gas_v_transf[:,0], gas_v_transf[:,1], gas_v_transf[:,2], gas_nf]))

	if i % 100 == 0:
		print(i)
	#plot spatial maps to see what's up 
		plt.scatter(star_r_transf[:,0], star_r_transf[:,1], alpha=.2, s=3)
		plt.scatter(0,0, c='r')
		plt.savefig('/Volumes/FRIEND/analogs/plots/TNG/rotation_maps/{}_star.png'.format(int(all_ids[i])))
		plt.close()
		plt.scatter(gas_r_transf[:,0], gas_r_transf[:,1], alpha=.2, s=3)
		plt.scatter(0,0, c='r')
		plt.savefig('/Volumes/FRIEND/analogs/plots/TNG/rotation_maps/{}_gas.png'.format(int(all_ids[i])))
		plt.close()



