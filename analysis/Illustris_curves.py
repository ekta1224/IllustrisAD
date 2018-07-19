import numpy as np 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 

'''
-do we try to pair the gas and stars?
-should we smooth the velocities?
'''

halo = '361428'

#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass_all, star_factor = np.loadtxt('../data/M31analog_{}_star_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True) 
gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass_all, gas_fraction = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True)

h = 0.704

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y, z):
	#convert to kpc (coordinates are already shifted to 0,0,0)
	x_kpc = x / h
	y_kpc = y / h
	z_kpc = z / h
	return np.sqrt((x_kpc**2) + (y_kpc**2) + (z_kpc**2)) #kpc

#calculate the radial velocity
def radial_v(x, y, z, vx, vy, vz):
	dot_product = (x * vx) + (y * vy) + (z * vz) 
	r_mag = distance(x, y, z)
	return dot_product / r_mag

#get rid of star partiles that are really wind
not_wind = star_factor > 0
star_x = star_x[not_wind]
star_y = star_y[not_wind]
star_z = star_z[not_wind]
star_vx = star_vx[not_wind]
star_vy = star_vy[not_wind]
star_vz = star_vz[not_wind]
star_factor = star_factor[not_wind]

star_vrad = radial_v(star_x, star_y, star_z, star_vx, star_vy, star_vz)

star_r = distance(star_x, star_y, star_z)
gas_r = distance(gas_x, gas_y, gas_z)

#only want particles/cells that correspond to a high HI gas fraction
neutral_gas = gas_fraction > 0.7 #what is the best limit?

gas_r = gas_r[neutral_gas]
gas_x = gas_x[neutral_gas]
gas_y = gas_y[neutral_gas]
gas_z = gas_z[neutral_gas]
gas_vx = gas_vx[neutral_gas]
gas_vy = gas_vy[neutral_gas]
gas_vz = gas_vz[neutral_gas]
gas_vrad = radial_v(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz)

#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
def vrot(x, y, z, vx, vy, vz):
	v_rad = radial_v(x, y, z, vx, vy, vz)
	v_tot = np.sqrt((vx**2) + (vy**2) + (vz**2))
	v_rot = np.sqrt((v_tot**2) + (v_rad**2))
	return v_rot

star_vrot = vrot(star_x, star_y, star_z, star_vx, star_vy, star_vz)
gas_vrot = vrot(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz)

#calculate the ages of stars from the scale factor
def star_age(scale_factor):
	zform = [(1./a) -1. for a in scale_factor]
	cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
	tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
	tage = np.array([13.8-t for t in tform])
	return tage

age = star_age(star_factor)

#divide the stars into 4 age bins
group1 = (.0 <= age) & (age <=.06) #average age 30 Myr
group2 = (.34 <= age) & (age <=.46) #average age 400 Myr
group3 = (1.7 <= age) & (age <=2.3) #average age 2 Gyr
group4 = (3.7 <= age) & (age <=4.3) #average age 4 Gyr

star1_r = star_r[group1]
star1_vrot = star_vrot[group1]

star2_r = star_r[group2]
star2_vrot = star_vrot[group2]

star3_r = star_r[group3]
star3_vrot = star_vrot[group3]

star4_r = star_r[group4]
star4_vrot = star_vrot[group4]

#divide stars and gas into radial and average the rotation velocities within each bin 
R_min = 0 #kpc
R_max = 21
delta_r= 0.25 #kpc

r_bins=np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)

#will contain the average rotation velocity at each radial bin
star1_avg_vrot = np.zeros(len(r_bins))
star2_avg_vrot = np.zeros(len(r_bins))
star3_avg_vrot = np.zeros(len(r_bins))
star4_avg_vrot = np.zeros(len(r_bins))
gas_avg_vrot = np.zeros(len(r_bins))

#dividing each age groups into radial bin
for i in range(len(r_bins)-1):
	star1_vrots=[b for b, a in zip(star1_vrot, star1_r) if a>=r_bins[i] and a<r_bins[i+1]]
	star1_avg_vrot[i]=np.median(star1_vrots)
	star2_vrots=[b for b, a in zip(star2_vrot, star2_r) if a>=r_bins[i] and a<r_bins[i+1]]
	star2_avg_vrot[i]=np.median(star2_vrots)
	star3_vrots=[b for b, a in zip(star3_vrot, star3_r) if a>=r_bins[i] and a<r_bins[i+1]]
	star3_avg_vrot[i]=np.median(star3_vrots)
	star4_vrots=[b for b, a in zip(star4_vrot, star4_r) if a>=r_bins[i] and a<r_bins[i+1]]
	star4_avg_vrot[i]=np.median(star4_vrots)
	gas_vrots=[b for b, a in zip(gas_vrot, gas_r) if a>=r_bins[i] and a<r_bins[i+1]]
	gas_avg_vrot[i]=np.median(gas_vrots)

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

star1_ad = va(gas_avg_vrot, star1_avg_vrot)
star2_ad = va(gas_avg_vrot, star2_avg_vrot)
star3_ad = va(gas_avg_vrot, star3_avg_vrot)
star4_ad = va(gas_avg_vrot, star4_avg_vrot)

#remove nan values (that resulted from having empty radial bins)
star1_ad = star1_ad[~np.isnan(star1_ad)]
star2_ad = star2_ad[~np.isnan(star2_ad)]
star3_ad = star3_ad[~np.isnan(star3_ad)]
star4_ad = star4_ad[~np.isnan(star4_ad)]

#plots-- 4 rotation curves and 1 histogram of AD
plt.scatter(r_bins[:-1], gas_avg_vrot[:-1], c = 'darkgrey', alpha = 0.4, label='gas')
plt.scatter(r_bins[:-1], star1_avg_vrot[:-1], c = 'b', alpha = 0.4, label='30 Myr')
plt.scatter(r_bins[:-1], star2_avg_vrot[:-1], c = 'm', alpha = 0.4, label='400 Myr')
plt.scatter(r_bins[:-1], star3_avg_vrot[:-1], c = 'k', alpha = 0.4, label='2 Gyr')
plt.scatter(r_bins[:-1], star4_avg_vrot[:-1], c = 'r', alpha = 0.4, label='4 Gyr')
plt.ylim(0,300)
plt.xlim(0,20)
plt.ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=13)
plt.xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
plt.legend(loc=4)
plt.savefig('{}_rotation_curve.pdf'.format(halo))
plt.close()

plt.hist(star1_ad[:-1], bins=range(-70, 70, 10), label='30 Myr', normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(star2_ad[:-1], bins=range(-70, 70, 10), label='400 Myr', normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m')
plt.hist(star3_ad[:-1], bins=range(-70, 70, 10), label='2 Gyr', normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='k')
plt.hist(star4_ad[:-1], bins=range(-70, 70, 10), label='4 Gyr', normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='r')
plt.legend(loc=1, frameon=False)
plt.xlim(-75,70)
plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
plt.savefig('{}_asymmetric_drift_hist.pdf'.format(halo))
