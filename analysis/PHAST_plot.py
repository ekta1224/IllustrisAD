'''
make Claire's age vs velocity dispersion plot for the PHAST proposal
'''

import numpy as np 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 
from astropy.coordinates import SkyCoord
from matplotlib import rc 

def single_plot():
	rc('font', family = 'serif')
	fig, ax=plt.subplots(1)
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()
	return

#read in the data for the analog
halo = 396140
star_x, star_y, star_z, star_vx, star_vy, star_vz, star_factor = np.loadtxt('/Volumes/Titan/analogs/IllustrisAD/data/M31analog_{}_star_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 3, 4, 5, 7), unpack = True) 

#Illustric value
h = 0.704

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y, z):
	#(coordinates are already shifted to 0,0,0)
	return np.sqrt((x**2) + (y**2) + (z**2)) #whatever units the input is in

#used for calculating velocities in the x y plane
def projected_distance(x, y, z):
	#(coordinates are already shifted to 0,0,0)
	return np.sqrt((x**2) + (y**2)) #whatever units the input is in

#calculate the radial velocity
def radial_v(x, y, z, vx, vy, vz): #just want in x y plane so removing z component
	#convert kpc/h to km
	x_km = x * 3.086e+16 / h
	y_km = y * 3.086e+16 / h
	z_km = z * 3.086e+16 / h
	dot_product = (x_km * vx) + (y_km * vy) 
	r_mag = projected_distance(x_km, y_km, z_km)
	return dot_product / r_mag #km/s

#get rid of star partiles that are really wind
not_wind = star_factor > 0
star_x = star_x[not_wind]
star_y = star_y[not_wind]
star_z = star_z[not_wind]
star_vx = star_vx[not_wind]
star_vy = star_vy[not_wind]
star_vz = star_vz[not_wind]
star_factor = star_factor[not_wind]

#get rid of particles that are too distant in the z direction
close = abs(star_z) / h <= 10#) & (star_y / h < 1 - 17 * (star_x / h) / 14) #second part is used for disk splitting tests
star_x = star_x[close]
star_y = star_y[close]
star_z = star_z[close]
star_vx = star_vx[close]
star_vy = star_vy[close]
star_vz = star_vz[close]
star_factor = star_factor[close]

#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
def vrot(x, y, z, vx, vy, vz): #just want planar so removing z component
	v_rad = radial_v(x, y, z, vx, vy, vz) #km/s
	v_tot = np.sqrt((vx**2) + (vy**2)) #km/s
	v_rot = np.sqrt((v_tot**2) - (v_rad**2))
	return v_rot

#smoothing ===============================================================
def smoothing(group, circle_size): #do not smooth the gas
	#group is going to encompass the age bin and the left or right side
	#get data for the correct age group
	star_xs = star_x[group]
	star_ys = star_y[group]
	star_zs = star_z[group]
	star_vxs = star_vx[group]
	star_vys = star_vy[group]
	star_vzs = star_vz[group]

	star_smoothed_vx = []
	star_smoothed_vy = []
	star_smoothed_vz = []
	star_goodcenter_x = []
	star_goodcenter_y = []
	star_goodcenter_z = []
	star_dispersion_Z = [] #vertical component
	star_dispersion_R = [] #radial component
	star_dispersion_phi = [] #azimuthal component
	c = SkyCoord(ra = star_xs / 13.86 / h, dec = star_ys / 13.86 / h, unit=(u.deg,u.deg))
	for i in range(len(star_xs)):
		c1 = SkyCoord(star_xs[i] / 13.86 / h, star_ys[i] / 13.86 / h, unit=(u.deg,u.deg)) 
		sep = c1.separation(c)
		good = sep.arcsecond <= circle_size #put stars into smoothing circle of this size
		star_x_velocities = star_vxs[good]
		star_y_velocities = star_vys[good]
		star_z_velocities = star_vzs[good]
		radial_vels = radial_v(star_xs[good], star_ys[good], star_zs[good], star_x_velocities, star_y_velocities, star_z_velocities) #to get dispersion
		vrots = vrot(star_xs[good], star_ys[good], star_zs[good], star_x_velocities, star_y_velocities, star_z_velocities) #to get dispersion)
		if len(star_x_velocities) >= 10:
			star_goodcenter_x.append(star_xs[i]) #kpc/h
			star_goodcenter_y.append(star_ys[i]) #kpc/h
			star_goodcenter_z.append(star_zs[i]) #kpc/h
			star_smoothed_vx.append(np.median(star_x_velocities)) #km/s
			star_smoothed_vy.append(np.median(star_y_velocities)) #km/s
			star_smoothed_vz.append(np.median(star_z_velocities)) #km/s
			star_dispersion_Z.append(np.std(star_z_velocities)) #km/s
			star_dispersion_R.append(np.std(radial_vels)) 
			star_dispersion_phi.append(np.std(vrots))

	return star_dispersion_Z #np.array((star_smoothed_vx)), np.array((star_smoothed_vy)), np.array((star_smoothed_vz)), np.array((star_goodcenter_x)), np.array((star_goodcenter_y)), np.array((star_goodcenter_z)), np.array((star_dispersion_Z)), np.array((star_dispersion_R)), np.array((star_dispersion_phi))

#calculate the ages of stars from the scale factor
def star_age(scale_factor):
	zform = [(1./a) -1. for a in scale_factor]
	cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
	tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
	tage = np.array([13.8-t for t in tform])
	return tage

age = star_age(star_factor)

#divide the stars into 4 age bins
group1 = (0.03 - 0.03 < age) & (age < 0.03 + 0.07) #Gyr; doing a dynamic range for the ages 0.0-0.06
group2 = (0.4 - 0.4 / 2 < age) & (age < 0.4 + 0.4 / 2) #0.2-0.6 
group3 = (2 - 2 / 2 < age) & (age < 2 + 2 / 2) #1-3
group4 = (4 - 4 / 4 < age) & (age < 4 + 4 / 4) #3-5

#divide the analog into two halves
left = star_y / h < -2.0093492768477392 * star_x / h
right = star_y / h > -2.0093492768477392 * star_x / h

#get smoothed data
age1_left_phis = smoothing(group1 * left, 275)
age1_right_phis = smoothing(group1 * right, 275)
age2_left_phis = smoothing(group2 * left, 200)
age2_right_phis = smoothing(group2 * right, 200)
age3_left_phis = smoothing(group3 * left, 200)
age3_right_phis = smoothing(group3 * right, 200)
age4_left_phis = smoothing(group4 * left, 200)
age4_right_phis = smoothing(group4 * right, 200)

#caluclate the median velocity dispersion
age1_left_phi = np.median(age1_left_phis)
age1_right_phi = np.median(age1_right_phis)
age2_left_phi = np.median(age2_left_phis)
age2_right_phi = np.median(age2_right_phis)
age3_left_phi = np.median(age3_left_phis)
age3_right_phi= np.median(age3_right_phis)
age4_left_phi = np.median(age4_left_phis)
age4_right_phi = np.median(age4_right_phis)

#plot it up
ages = [0.03, 0.4, 2, 4] #Gyrs
lefts = [age1_left_phi, age2_left_phi, age3_left_phi, age4_left_phi]
rights = [age1_right_phi, age2_right_phi, age3_right_phi, age4_right_phi]

single_plot()
#plt.plot(ages, lefts, 'bo', label='left half')
#plt.plot(ages, rights, 'ro', label='right half')
plt.scatter(ages[0], lefts[0], marker='s', c='b')
plt.scatter(ages[1], lefts[1], marker='s', c='k', label='left half')
plt.scatter(ages[2], lefts[2], marker='s', c='purple')
plt.scatter(ages[3], lefts[3], marker='s', c='r')
plt.scatter(ages[0], rights[0], marker='o', c='b')
plt.scatter(ages[1], rights[1], marker='o', c='k', label='right half')
plt.scatter(ages[2], rights[2], marker='o', c='purple')
plt.scatter(ages[3], rights[3], marker='o', c='r')
plt.xlabel(r'$\rm Mean\ Age\ (Gyr)$', fontsize=13)
plt.ylabel(r'$\rm Median\ \sigma_{v}\ (km/s)$', fontsize=13)
plt.legend(frameon=False)
plt.savefig('/Users/amandaquirk/Desktop/PHAST_vel_disp.pdf', bbox_inches='tight')
plt.close()

# single_plot()
# plt.hist(age2_left_phis, bins=range(0, 100, 5), label='left, #=' + '{}'.format(len(age2_left_phis)), normed=1, color='b')
# plt.hist(age2_right_phis, bins=range(0, 100, 5), label='right, #=' + '{}'.format(len(age2_right_phis)), normed=1, color='r', alpha=0.35)
# plt.legend(frameon=False)
# plt.savefig('/Users/amandaquirk/Desktop/age_2_vel_disp.png')
# plt.close()






