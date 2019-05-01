import numpy as np 
import astropy.units as u 
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt 

'''
purpose: find analogs that do not have enough neutral gas to be reliable in our AD study
'''

halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt')

# #first I look at the mass and number of particles of the NEUTRAL gas ====================================================
# neutral_fraction = [.7, .75, .8, .85, .9, .95]

# neutral_mass_7 = np.zeros_like(halos)
# num_neutral_particles_7 = np.zeros_like(halos)
# neutral_mass_75 = np.zeros_like(halos)
# num_neutral_particles_75 = np.zeros_like(halos)
# neutral_mass_8 = np.zeros_like(halos)
# num_neutral_particles_8 = np.zeros_like(halos)
# neutral_mass_85 = np.zeros_like(halos)
# num_neutral_particles_85 = np.zeros_like(halos)
# neutral_mass_9 = np.zeros_like(halos)
# num_neutral_particles_9 = np.zeros_like(halos)
# neutral_mass_95 = np.zeros_like(halos)
# num_neutral_particles_95 = np.zeros_like(halos)

# for i in range(len(halos)):
# 	print(halos[i])
# 	#coordinates ckpc/h -- since at z=0, kpc/h
# 	#speeds km sqrt(a)/s -- since at z=0, km/s
	
# 	#Illustris value
# 	h = 0.704

# 	def gas_smoothing(halo, fraction):
# 		gas_x, gas_y, gas_z, gas_mass, gas_fraction = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 6, 7,), unpack = True)
# 		#only want particles/cells that correspond to a high HI gas fraction
# 		close_neutral_gas = (gas_fraction >= fraction) & (gas_z / h <=10) #what is the best limit?
	
# 		gas_x = gas_x[close_neutral_gas]
# 		gas_y = gas_y[close_neutral_gas]
# 		gas_z = gas_z[close_neutral_gas]
# 		gas_mass = gas_mass[close_neutral_gas]
		
# 		#smoothing gas (for tests)
# 		gas_smoothed_mass = []
# 		c = SkyCoord(ra = gas_x / 13.86 / h, dec = gas_y / 13.86 / h, unit=(u.deg,u.deg))
# 		for i in range(len(gas_x)):
# 			c1 = SkyCoord(gas_x[i] / 13.86 / h, gas_y[i] / 13.86 / h, unit=(u.deg,u.deg)) 
# 			sep = c1.separation(c)
# 			good = sep.arcsecond <= 275 #put stars into smoothing circle of this size
# 			if sum(good) >= 10:
# 				gas_smoothed_mass.append(gas_mass[i] * 10**10 / h) #Msun
# 		return sum(gas_smoothed_mass), len(gas_smoothed_mass)

# 	neutral_mass_7[i], num_neutral_particles_7[i] = gas_smoothing(halos[i], .7)
# 	print("done with .7")
# 	neutral_mass_75[i], num_neutral_particles_75[i] = gas_smoothing(halos[i], .75)
# 	print("done with .75")
# 	neutral_mass_8[i], num_neutral_particles_8[i] = gas_smoothing(halos[i], .8)
# 	print("done with .8")
# 	neutral_mass_85[i], num_neutral_particles_85[i] = gas_smoothing(halos[i], .85)
# 	print("done with .85")
# 	neutral_mass_9[i], num_neutral_particles_9[i] = gas_smoothing(halos[i], .9)
# 	print("done with .9")
# 	neutral_mass_95[i], num_neutral_particles_95[i] = gas_smoothing(halos[i], .95)

# np.savetxt('/Volumes/FRIEND/analogs/data/neutral_gas_mass_all.txt', np.c_[halos, neutral_mass_7, num_neutral_particles_7, neutral_mass_75, num_neutral_particles_75, neutral_mass_8, num_neutral_particles_8, neutral_mass_85, num_neutral_particles_85, neutral_mass_9, num_neutral_particles_9, neutral_mass_95, num_neutral_particles_95], fmt='%1.16f', delimiter=' ', header='ID, total neutral gas mass and number of neutral gas particles for: neutral fraction= .7, .75, .8, .85, .9, .95')

# below looks at the fraction of neutral gas to the total gas =======================================================

#neutral_gas_mass7, neutral_gas_mass75, neutral_gas_mass8, neutral_gas_mass85, neutral_gas_mass9, neutral_gas_mass95 = np.loadtxt('/Volumes/FRIEND/analogs/data/gas_mass_all.txt', usecols=(1, 3, 5, 7, 9, 11), unpack=True)
total_mass = np.zeros_like(halos)

for i in range(len(halos)):
	#print(halos[i])
	gas_x, gas_y, gas_z, gas_mass = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(int(halos[i])), usecols=(0, 1, 2, 6,), unpack = True)
	
	#coordinates ckpc/h -- since at z=0, kpc/h
	#speeds km sqrt(a)/s -- since at z=0, km/s
	
	#Illustris value
	h = 0.704

	#only want particles/cells that correspond to a high HI gas fraction
	close_gas = gas_z / h <=10 #what is the best limit?
	
	gas_mass = gas_mass[close_gas] * 10**10 / h
	total_mass[i] = (sum(gas_mass))

def fraction_hist(mass_array, fraction):
	plt.hist(mass_array/total_mass, bins=np.linspace(-0.01, 0.6, 20))
	plt.xlabel('Neutral Gas Mass/ Total Gas Mass')
	plt.savefig('/Users/amandaquirk/Desktop/{}_fraction.png'.format(fraction))
	plt.close()
	return

neutral_mass_7, neutral_mass_75, neutral_mass_8, neutral_mass_85, neutral_mass_9, neutral_mass_95 = np.loadtxt('/Volumes/FRIEND/analogs/data/neutral_gas_mass_all.txt', usecols=(1, 3, 5, 7, 9, 11), unpack=True)

# fraction_hist(neutral_mass_7, 7)
# fraction_hist(neutral_mass_75, 75)
# fraction_hist(neutral_mass_8, 8)
# fraction_hist(neutral_mass_85, 85)
# fraction_hist(neutral_mass_9, 9)
# fraction_hist(neutral_mass_95, 95)

print('Fraction = 0.7==============================')
fraction_7 = neutral_mass_7 / total_mass
for i in range(len(fraction_7)):
	if fraction_7[i] < 0.25:
		print(halos[i])	

print('Fraction = 0.75==============================')
fraction_75 = neutral_mass_75 / total_mass
for i in range(len(fraction_7)):
	if fraction_75[i] < 0.25:
		print(halos[i])	

print('Fraction = 0.8==============================')
fraction_8 = neutral_mass_8 / total_mass
for i in range(len(fraction_7)):
	if fraction_8[i] < 0.25:
		print(halos[i])	

print('Fraction = 0.85==============================')
fraction_85 = neutral_mass_85 / total_mass
for i in range(len(fraction_7)):
	if fraction_85[i] < 0.25:
		print(halos[i])	

print('Fraction = 0.9==============================')
fraction_9 = neutral_mass_9 / total_mass
for i in range(len(fraction_7)):
	if fraction_9[i] < 0.25:
		print(halos[i])	

print('Fraction = 0.95==============================')
fraction_95 = neutral_mass_95 / total_mass
for i in range(len(fraction_7)):
	if fraction_95[i] < 0.25:
		print(halos[i])	






