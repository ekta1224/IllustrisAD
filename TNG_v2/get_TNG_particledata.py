import requests 
import numpy as np 

#function that gives path to data
def get(path, params=None):
	# make HTTP GET request to path
	headers = {"api-key":"673159e0bdccc360baf42bfb32017ec7"}
	r = requests.get(path, params=params, headers=headers)
	
	# raise exception if response code is not HTTP SUCCESS (200)
	r.raise_for_status()
	
	if r.headers['content-type'] == 'application/json':
	    return r.json() # parse json responses automatically
	
	if 'content-disposition' in r.headers:
	    filename = r.headers['content-disposition'].split("filename=")[1]
	    with open(filename, 'wb') as f:
	        f.write(r.content)
	    return filename # return the filename string
	
	return r

# #step 1 -- do a mass cut ======================================================================
# mass_min = 1 * 10**12 / 1e10 * 0.704                                                            
# mass_max = 2 * 10**12 / 1e10 * 0.704                                                           

# search_query = "?mass__gt=" + str(mass_min) + "&mass__lt=" + str(mass_max)                     
# url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/" + search_query      
# subhalos = get(url)
# all_subhalos = get(url, {'limit':subhalos['count']}) 
# print('{} analogs passed the initial mass cut'.format(subhalos['count']))

# #step 2 -- only want halos that are the center of their FoF group ============================
# #central_flag = [] #0 = npt center and 1 = center
# ids = [all_subhalos['results'][i]['id'] for i in range(all_subhalos['count'])] 
# centers_id = [] #contains the IDs of the halos that are the center of ther FoF group
# for i in range(len(ids)):#all_subhalos['count']):
# 	if i % 100 == 0:
# 		print(i)
# 	sub = get(all_subhalos['results'][i]['url'])
# 	#central_flag.append(sub['primary_flag'])
# 	if sub['primary_flag'] == 1:
# 		centers_id.append(ids[i])
# print('{} analogs passed the FoF center cut'.format(len(centers_id)))
# np.savetxt('/Users/amandaquirk/Desktop/TNG_central_IDs.txt', centers_id)

# #step 3 -- look at the satellites around the centrals ==========================================
# #read in IDs of central FoF massive halos, since step 2 isn't needed every time 
# ids = np.loadtxt('/Users/amandaquirk/Desktop/TNG_central_IDs.txt') #ids of the massive centers

# #will contain the ids of the subhalos that either are isolated or are have 1 M33-like companion
# isolated = []
# M33_comp = []

# h = 0.704
# def distance(centerx, centery, centerz, satx, saty, satz): #code units of ckpc / h
# 	dist = np.sqrt((centerx - satx)**2 + (centery - saty)**2 + (centerz - satz)**2)
# 	return dist / h #kpc

# for i in range(len(ids)): 
# 	if i % 100 == 0:
#  		print(i)
# 	#get info about the central 
# 	url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/" + str(int(ids[i]))
# 	subhalo = get(url)
# 	group_num = subhalo['grnr']
# 	central_x = subhalo['pos_x']
# 	central_y = subhalo['pos_y']
# 	central_z = subhalo['pos_z']

# 	#get all of the subhalos in the central's friend of friend group
# 	fof_url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/?grnr=" + str(group_num)
# 	fof_group = get(fof_url)
# 	fof_group = get(fof_url, {'limit':fof_group['count']})
# 	masses = [fof_group['results'][j]['mass_log_msun'] for j in range(fof_group['count'])] #log(Msun) ordered largest to smallest
# 	central_mass = 10**masses[0] / h #Msun
# 	rvir = 260 * (central_mass / 10**12)**(1 / 3) #kpc, following patel et al 2017a

# 	if len(masses) == 1:
# 		isolated.append(ids[i])

# 	else:
# 		#get position of largest satellite
# 		id_large_sat = fof_group['results'][1]['id'] #id of the largest satellite
# 		large_sat_url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/" + str(id_large_sat)
# 		large_sat = get(large_sat_url)
# 		large_sat_x = large_sat['pos_x']
# 		large_sat_y = large_sat['pos_y']
# 		large_sat_z = large_sat['pos_z']
# 		large_sat_dist = distance(central_x, central_y, central_z, large_sat_x, large_sat_y, large_sat_z) #kpc
	
# 		#find analogs that fit our criteria
# 		if masses[1] < np.log10(8e10): #WHEN BACK FROM TURKEY NEED TO CONFIRM THIS CRITERIA WITH EKTA
# 			isolated.append(ids[i])
# 		elif (masses[1] / h < np.log10(32e10)) & (masses[1] / h > np.log10(8e10)) & (masses[2] / h < np.log10(8e10)) & (large_sat_dist < rvir):
# 			M33_comp.append(ids[i])

# print('We started with {} central halos'.format(len(ids)))
# print('{} of these are isolated'.format(len(isolated)))
# print('{} of these have M33-like companions'.format(len(M33_comp)))
# np.savetxt('/Users/amandaquirk/Desktop/isolated_central_IDs.txt', isolated)
# np.savetxt('/Users/amandaquirk/Desktop/M33com_central_IDs.txt', M33_comp)

#step 4 -- grab the halo and particle data for all of the analogs found in the steps above 
import h5py 

#read in IDs
#noM33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_noM33_TNG100.txt')
#M33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_M33_TNG100.txt')
revised_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31analogs_halo_props_TNG100_revised.txt', usecols=(0,), unpack=True)
revised_ids = [int(a) for a in revised_ids]

#list the parameters I want
def save_particle_data(subhalo_id):
	params = {'stars':'Coordinates,Velocities,GFM_StellarFormationTime,Masses,GFM_Metallicity','gas':'Coordinates,Velocities,Masses,NeutralHydrogenAbundance,StarFormationRate,GFM_Metallicity','dm':'Coordinates,Velocities,Potential'} 
	url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/" + str(int(subhalo_id))
	saved_filename = get(url + "/cutout.hdf5",params) #hdf5 file with all the data
	f = h5py.File(saved_filename, 'r') #read the file in order to save the data
	# stars = f['PartType4']
	# gas = f['PartType0']
	# dm = f['PartType1']
	star_x = f['PartType4']['Coordinates'][:,0]
	star_y = f['PartType4']['Coordinates'][:,1]
	star_z = f['PartType4']['Coordinates'][:,2]
	star_vx = f['PartType4']['Velocities'][:,0]
	star_vy = f['PartType4']['Velocities'][:,1]
	star_vz = f['PartType4']['Velocities'][:,2]
	star_mass = f['PartType4']['Masses']
	star_tform = f['PartType4']['GFM_StellarFormationTime']
	star_met = f['PartType4']['GFM_Metallicity']
	
	gas_x = f['PartType0']['Coordinates'][:,0]
	gas_y = f['PartType0']['Coordinates'][:,1]
	gas_z = f['PartType0']['Coordinates'][:,2]
	gas_vx = f['PartType0']['Velocities'][:,0]
	gas_vy = f['PartType0']['Velocities'][:,1]
	gas_vz = f['PartType0']['Velocities'][:,2]
	gas_mass = f['PartType0']['Masses']
	gas_nf = f['PartType0']['NeutralHydrogenAbundance']
	gas_sfr = f['PartType0']['StarFormationRate']
	gas_met = f['PartType0']['GFM_Metallicity']
	
	dm_x = f['PartType1']['Coordinates'][:,0]
	dm_y = f['PartType1']['Coordinates'][:,1]
	dm_z = f['PartType1']['Coordinates'][:,2]
	dm_vx = f['PartType1']['Velocities'][:,0]
	dm_vy = f['PartType1']['Velocities'][:,1]
	dm_vz = f['PartType1']['Velocities'][:,2]
	dm_mass = np.ones(len(dm_vz)) * 6.3 * 10**6
	
	#save the data to text files -- sorry Bruno 
	np.savetxt('/Volumes/Titan/analogs/TNGdata/particledata/{}_star_properties_TNGv2.txt'.format(subhalo_id), np.c_[star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass, star_tform, star_met])
	np.savetxt('/Volumes/Titan/analogs/TNGdata/particledata/{}_gas_properties_TNGv2.txt'.format(subhalo_id), np.c_[gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass, gas_nf, gas_sfr, gas_met])
	np.savetxt('/Volumes/Titan/analogs/TNGdata/particledata/{}_dm_properties_TNGv2.txt'.format(subhalo_id), np.c_[dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz, dm_mass])
	return

#for subhalo in revised_ids:
#	save_particle_data(subhalo)
save_particle_data(int(4.863410000000000000e+05))

# def get_subhalo_params(subhalo_id): #can I add merger props here??
# 	url = "http://www.tng-project.org/api/TNG100-1/snapshots/z=0/subhalos/" + str(int(subhalo_id))
# 	subhalo = get(url)
# 	x = subhalo['pos_x']
# 	y = subhalo['pos_y']
# 	z = subhalo['pos_z']
# 	vx = subhalo['vel_x']
# 	vy = subhalo['vel_y']
# 	vz = subhalo['vel_z']
# 	mass = subhalo['mass']
# 	radius = subhalo['halfmassrad']
# 	stellar_mass = subhalo['mass_stars']

# 	return x, y, z, vx, vy, vz, mass, radius, stellar_mass 

# all_ids = list(M33_ids) + list(noM33_ids)
# # all_x = []
# # all_y = []
# # all_z = []
# # all_vx = []
# # all_vy = []
# # all_vz = []
# # all_mass = []
# # all_radii = []
# #all_smass = []
# for id in noM33_ids:
# 	save_particle_data(id)
# 	data = get_subhalo_params(id)
# 	all_x.append(data[0])
# 	all_y.append(data[1])
# 	all_z.append(data[2])
# 	all_vx.append(data[3])
# 	all_vy.append(data[4])
# 	all_vz.append(data[5])
# 	all_mass.append(data[6])
# 	all_radii.append(data[7])
# 	all_smass.append(data[8])

# np.savetxt('/Volumes/FRIEND/analogs/TNGdata/TNG_halo_props.txt', np.c_[all_ids, all_x, all_y, all_z, all_vx, all_vy, all_vz, all_mass, all_radii, all_smass])






















