import numpy as np 
import matplotlib.pyplot as plt 

'''
1. CHECK UNITS FOR EVERYTHING-- position and velocity
2. do first without any smoothing but smoothing might be needed (maybe the radial binning will remove enough scatter)
'''

#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
star_x, star_y, star_z, star_v, star_mass_all, star_age = np.loadtxt('star_file.txt', unpack = True) 
gas_x, gas_y, gas_z, gas_v, gas_mass_all, gas_fraction = np.loadtxt('gas_file.txt', unpack = True)
dm_x, dm_y, dm_z, dm_mass_all = np.loadtxt('dm_file.txt', unpack = True)

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y, z):
	return np.sqrt((x**2) + (y**2) + (z**2))

star_r_all = distance(stars_x, stars_y, stars_z)
gas_r_all = distance(gas_x, gas_y, gas_z)
dm_r_all = distance(dm_x, dm_y, dm_z)

#combines the mass and radii for all particles-- to be used in the vrot function
mass = star_mass_all + gas_mass_all + dm_mass_all #make sure this combines the lists and doesn't add them
radii = star_r_all + gas_r_all + dm_r_all #make sure this combines the lists and doesn't add them

#only want particles/cells that correspond to a high HI gas fraction
neutral_gas = 0.7

star_r = star_r_all[neutral_gas]
star_v = star_v[neutral_gas]
star_age = star_age[neutral_gas]
gas_r = gas_r_all[neutral_gas]
gas_v = gas_v[neutral_gas]

#calculate the ages of stars
age = #Ekta's scale factor to age function

#divide the stars into 4 age bins
group1 = (.026 <= age) & (age <=.034) #average age 30 Myr
group2 = (.36 <= age) & (age <=.44) #average age 400 Myr
group3 = (1.6 <= age) & (age <=2.4) #average age 2 Gyr
group4 = (3.6 <= age) & (age <=4.4) #average age 4 Gyr

star1_r = star_r[group1]
star1_v = star_v[group1]

star2_r = star_r[group2]
star2_v = star_v[group2]

star3_r = star_r[group3]
star3_v = star_v[group3]

star4_r = star_r[group4]
star4_v = star_v[group4]

#calculate the enclosed mass
def enclosed_mass(r, all_masses, all_radii):
	masses = [a for a,b in zip(all_star_masses, all_star_radii) if b < r]
	return np.sum(masses)

#calculate the rotation speed, assuming spherical
def vrot(r, all_masses, all_radii):
	M = enclosed_mass(r, all_masses, all_radii)
	G = 4.302* 10**(-3) #pc M_sun^-1 (km/s)^2
	return np.sqrt(G * M / r)

#divide stars and gas into radial bins and average the velocities within each bin 

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

#plots-- 4 rotation curves and 1 histogram of AD




