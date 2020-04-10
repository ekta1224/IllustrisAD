#import requests 
import numpy as np
import illustris_python as il
#import scipy.linalg as la
#from correct_position import *

#function that gives path to data
# def get(path, params=None):
# 	# make HTTP GET request to path
# 	headers = {"api-key":"673159e0bdccc360baf42bfb32017ec7"}
# 	r = requests.get(path, params=params, headers=headers)
	
# 	# raise exception if response code is not HTTP SUCCESS (200)
# 	r.raise_for_status()
	
# 	if r.headers['content-type'] == 'application/json':
# 	    return r.json() # parse json responses automatically
	
# 	if 'content-disposition' in r.headers:
# 	    filename = r.headers['content-disposition'].split("filename=")[1]
# 	    with open(filename, 'wb') as f:
# 	        f.write(r.content)
# 	    return filename # return the filename string
	
# 	return r

# not doing the commented out below lines since downloaded the data
# snap_url = "http://www.tng-project.org/api/TNG100-1/snapshots/99/"
# snap_info = get(snap_url)
# num_groups = range(snap_info['num_groups_fof']) #numbers of FoF groups in z=0

# ================================================================================================================================
#step 1 -- grab the raw data I downloaded UGH can you believe this is actually the faster way
basePath = '/Volumes/Titan/analogs/TNGdata/group_cat'
halos = il.groupcat.loadHalos(basePath,99,fields=['GroupFirstSub', 'GroupPos', 'GroupVel', 'Group_M_TopHat200', 'Group_R_TopHat200']) #positions, velocities, and radius all in comoving/code units
print('Read in the halos')
subhalos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMassType', 'SubhaloVmax', 'SubhaloPos', 'SubhaloVel']) #mass in code units, vmax in km/s
print('Read in the subhalos')

#step 2 -- find halos with correct halo mass
h = 0.704 #0.72
mass_min = 1 * (10**12 / 1e10) * h #7.30 * 10*11                                                            
mass_max = 2 * (10**12 / 1e10) * h #2.21 * 10*12
halo_mass = halos['Group_M_TopHat200']

halo_cut = (mass_min < halo_mass) & (halo_mass < mass_max) 
# analog_pos = halos['GroupPos'][halo_cut]
# analog_vel = halos['GroupVel'][halo_cut]
analog_halo_mass = halo_mass[halo_cut]
#analog_r = halos['Group_R_TopHat200'][halo_cut]
first_sub = halos['GroupFirstSub'][halo_cut]

#grab some data for these analogs from the subhalo catalogue
sub_stellar_mass = subhalos['SubhaloMassType'].T[4] #stars are particle 4
sub_stellar_mass = sub_stellar_mass[first_sub]
analog_pos = subhalos['SubhaloPos'][first_sub]
analog_vel = subhalos['SubhaloVel'][first_sub]
sub_vmax = subhalos['SubhaloVmax'][first_sub] #km/s, no weird unit
print('Did the halo mass cut = {}'.format(sum(halo_cut)))

#step 3 -- look at subhalo catalogue and find analogs that pass a stellar mass and vmax cut
M31_mass_min = 5 * (10**10 / 1e10) * h                                                            
M31_mass_max = 2 * (10**11 / 1e10) * h 
Vmax_lim = 0 #km/s

M31_analog = (M31_mass_min < sub_stellar_mass) & (sub_stellar_mass < M31_mass_max) & (Vmax_lim < sub_vmax)
analog_pos = analog_pos[M31_analog]
analog_vel = analog_vel[M31_analog]
analog_halo_mass = analog_halo_mass[M31_analog]
halo_ID = first_sub[M31_analog]
sub_stellar_mass = sub_stellar_mass[M31_analog]
sub_vmax = sub_vmax[M31_analog]
print('Did the stellar mass and vmax cut = {}'.format(sum(M31_analog)))

#step 4 -- save data of the halos
np.savetxt('/Users/amandaquirk/Desktop/M31analogs_halo_props_TNG100_revised.txt', np.c_[halo_ID, analog_pos[:,0], analog_pos[:,1], analog_pos[:,2], analog_vel[:,0], analog_vel[:,1], analog_vel[:,2], analog_halo_mass, sub_stellar_mass, sub_vmax], header='ID, analog position x y z (ckpc/h), analog velocity vx vy vz (km/s/a), halo mass (10^10 Msun/h), stellar mass (10^10 Msun/h), vmax (km/s)')  

'''
==================================================================================================================================
below was for the original analogs -- not doing M33 component for reviewer's comments
'''
# #step 2 -- look at the satellites of each M31-like halo 
# M33 = []
# no_M33 = []
# history = 0
# no_history = 0
# print('Checking subhalos of each M31 analog')
# for i in range(len(first_sub)):
# 	#if i % 10**2 == 0:
# 		#print(i)
# 	mvir = analog_mass[i]
# 	rvir = analog_r[i]
# 	first_sub_ind = first_sub[i]
# 	#go through the 50 most massive subhalos of each FoF group
# 	sub_ids = []
# 	sub_past_masses = [] 
# 	#print('fs={}'.format(first_sub_ind))
# 	for j in np.arange(1, 50, 1): #only want to look at the first 50; this is already more sufficient than we need to be probably -- note to self, the SubField catalog is organized by subhalo ID; the subhalos are organized roughly by FoF group,so by looking at the subhalos directly after the First Sub in the catalog, we're probing ones that could be a part of that halo (ensuring the subhalo is within the Rvir ensures we don't include halos in another FoF group)
# 		if la.norm(correct_position(sub_pos[first_sub_ind + j], sub_pos[first_sub_ind])) < rvir: #subhalo within Rvir of halo
# 			try:
# 				#print('j={}'.format(j))
# 				branch_masses = il.sublink_gal.loadTree(basePath, 99, first_sub_ind + j, fields=['SubhaloMass'], onlyMPB=False) #pull the subhalo's merger tree and grab its past masses
# 				mass_max = np.max(branch_masses)
# 				#print('mm={}'.format(mass_max))
# 				sub_ids.append(first_sub_ind + j) #this subhalo's ID
# 				if mass_max is not None: #not all missing histories flagged as errors (warnings instead) 
# 					sub_past_masses.append(mass_max) #the maximum mass of one particular subhalo
# 					#sub_ids.append(first_sub_ind + j) #this subhalo's ID
# 					history += 1
# 				else:
# 					no_history += 1
# 					sub_past_masses.append(sub_mass[first_sub_ind + j])
# 			except AttributeError:
# 				print('No merger history, skip and keep going') #if no merger history, can't find max mass
# 				continue
	
# 	#now that we have the max mass of all(ish) the subhalos in a halo, we want to grab the mass and ID of the most massive one
# 	sorted_ids = np.array(sub_ids)[np.argsort(sub_past_masses)]
# 	mass_sat_id = sorted_ids[-1]
# 	mass_sat_m = np.max(sub_past_masses)
# 	current_sat_m = sub_mass[mass_sat_id] #z=0 mass
# 	#check is a halo has a M33-like companion (max mass matches M33 estimate AND z=0 mass clears a minimum)
# 	if (8 < mass_sat_m < 32) and (current_sat_m / h > 1):
# 		M33.append(first_sub_ind)
# 	else:
# 		no_M33.append(first_sub_ind)

# #save two files: one with lists of analogs with M33 companions and one with a list of analogs that do not have a companion
# np.savetxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_M33_TNG100.txt', M33)
# np.savetxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_noM33_TNG100.txt', no_M33)
# #print(history, no_history)