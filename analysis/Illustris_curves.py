import numpy as np 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 
from astropy.coordinates import SkyCoord
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

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

'''
-divides stellar particles into 4 age bins
-eliminates gas particles that don't have a 'high' fraction of neutral hydrogen
-smooths partcile velocities to better match the observation analysis and to create spatial maps
-construct a rotation curve for the gas and star particles, binned radially
-calculates the asymmetric drift at each radial pin and outputs a histogram to compare the four age bins
'''

halos = np.loadtxt('/Volumes/Titan/analogs/TNGdata/M31analogs_halo_props_TNG100_revised.txt', usecols=(0,), unpack=True)
halos = [int(a) for a in halos]
#halos = [419510]# [3.564270000000000000e+05, 3.874560000000000000e+05, 4.000040000000000000e+05, 4.053000000000000000e+05, 4.089160000000000000e+05, 4.122090000000000000e+05, 4.174200000000000000e+05, 4.236530000000000000e+05, 4.282060000000000000e+05, 4.341020000000000000e+05]
# noM33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_noM33_TNG100.txt')
# M33_ids = np.loadtxt('/Users/amandaquirk/Desktop/M31_analogs_IDs_M33_TNG100.txt')
# all_ids = list(noM33_ids) #+ list(M33_ids)
#halos = [int(a) for a in halos]

#below contains the median of the velocity dispersions for each component for each group
star1_med_disp_Z = []
star2_med_disp_Z = []
star3_med_disp_Z = []
star4_med_disp_Z = []
star1_med_disp_R = []
star2_med_disp_R = []
star3_med_disp_R = []
star4_med_disp_R = []
star1_med_disp_phi = [] 
star2_med_disp_phi = [] 
star3_med_disp_phi = [] 
star4_med_disp_phi = [] 

for halo in halos:
	print(halo)
	#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
	#AMANDA CHECK THE COLUMNS FOR AGE AND NEUTRAL FRACTION; DIFFERENT FOR TNG AND ILLUSTRIS FILES
	star_x, star_y, star_z, star_vx, star_vy, star_vz, star_factor = np.loadtxt('/Volumes/Titan/analogs/TNGdata/rotated//{}_star_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 3, 4, 5, 6), unpack = True) #/Volumes/FRIEND/analogs/TNGdata/{}_rotated_star_particles.tx
	gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_fraction = np.loadtxt('/Volumes/Titan/analogs/TNGdata/rotated/{}_gas_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 3, 4, 5, 6), unpack = True)
	
	#coordinates ckpc/h -- since at z=0, kpc/h
	#speeds km sqrt(a)/s -- since at z=0, km/s
	
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
	
	star_vrad_all = radial_v(star_x, star_y, star_z, star_vx, star_vy, star_vz) #km/s
	star_r_all = distance(star_x / h, star_y / h, star_z / h) #kpc
	
	#only want particles/cells that correspond to a high HI gas fraction
	# plt.clf()
	# plt.hist(gas_fraction)
	# plt.savefig('/Users/amandaquirk/Desktop/HI_fraction_{}.png'.format(int(halo)))
	# plt.close()
	neutral_gas = gas_fraction > 0.75 #what is the best limit?
	
	gas_x = gas_x[neutral_gas]
	gas_y = gas_y[neutral_gas]
	gas_z = gas_z[neutral_gas]
	gas_vx = gas_vx[neutral_gas]
	gas_vy = gas_vy[neutral_gas]
	gas_vz = gas_vz[neutral_gas]

	#get rid of particles that are too distant in the z direction
	close = abs(gas_z) / h <= 10 #) & (gas_y / h < 1 - 17 * (gas_x / h) / 14) #second part is used for disk splitting tests
	gas_x = gas_x[close]
	gas_y = gas_y[close]
	gas_z = gas_z[close]
	gas_vx = gas_vx[close]
	gas_vy = gas_vy[close]
	gas_vz = gas_vz[close]
	
	gas_vrad = radial_v(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz) #km/s
	gas_r = distance(gas_x / h, gas_y / h, gas_z / h) #kpc

	#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
	def vrot(x, y, z, vx, vy, vz): #just want planar so removing z component
		v_rad = radial_v(x, y, z, vx, vy, vz) #km/s
		v_tot = np.sqrt((vx**2) + (vy**2)) #km/s
		v_rot = np.sqrt((v_tot**2) - (v_rad**2))
		return v_rot
	
	#smoothing ===============================================================
	
	def smoothing(group, circle_size): #do not smooth the gas
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

		return np.array((star_smoothed_vx)), np.array((star_smoothed_vy)), np.array((star_smoothed_vz)), np.array((star_goodcenter_x)), np.array((star_goodcenter_y)), np.array((star_goodcenter_z)), np.array((star_dispersion_Z)), np.array((star_dispersion_R)), np.array((star_dispersion_phi))

	def gas_smoothing(circle_size):
		#smoothing gas (for tests)
		gas_smoothed_vx = []
		gas_smoothed_vy = []
		gas_smoothed_vz = []
		gas_goodcenter_x = []
		gas_goodcenter_y = []
		gas_goodcenter_z = []
		gas_dispersion = []
		c = SkyCoord(ra = gas_x / 13.86 / h, dec = gas_y / 13.86 / h, unit=(u.deg,u.deg))
		for i in range(len(gas_x)):
			c1 = SkyCoord(gas_x[i] / 13.86 / h, gas_y[i] / 13.86 / h, unit=(u.deg,u.deg)) 
			sep = c1.separation(c)
			good = sep.arcsecond <= circle_size #put stars into smoothing circle of this size
			gas_x_velocities = gas_vx[good]
			gas_y_velocities = gas_vy[good]
			gas_z_velocities = gas_vz[good]
			radial_vels = radial_v(gas_x[good], gas_y[good], gas_z[good], gas_x_velocities, gas_y_velocities, gas_z_velocities) #to get dispersion
			if len(gas_x_velocities) >= 10:
				gas_goodcenter_x.append(gas_x[i]) #kpc/h
				gas_goodcenter_y.append(gas_y[i]) #kpc/h
				gas_goodcenter_z.append(gas_z[i]) #kpc/h
				gas_smoothed_vx.append(np.median(gas_x_velocities)) #km/s
				gas_smoothed_vy.append(np.median(gas_y_velocities)) #km/s
				gas_smoothed_vz.append(np.median(gas_z_velocities)) #km/s
				gas_dispersion.append(np.std(gas_z_velocities))
	
		return np.array((gas_smoothed_vx)), np.array((gas_smoothed_vy)), np.array((gas_smoothed_vz)), np.array((gas_goodcenter_x)), np.array((gas_goodcenter_y)), np.array((gas_goodcenter_z)), np.array((gas_dispersion))
	
	#=========================================================================
	'''
	keep if not smoothing
	star_vrot = vrot(star_x, star_y, star_z, star_vx, star_vy, star_vz)
	'''
	
	gas_vrot = vrot(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz)
	
	#calculate the ages of stars from the scale factor
	def star_age(scale_factor):
		zform = [(1./a) -1. for a in scale_factor]
		cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
		tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
		tage = np.array([13.8-t for t in tform])
		return tage
	
	age = star_age(star_factor)
	# plt.hist(age)
	# plt.savefig('/Users/amandaquirk/Desktop/{}_stellar_age.png'.format(halo))
	# plt.close()
	
	#divide the stars into 4 age bins
	group1 = age <= 1 #Gyr
	group2 = (1 < age) & (age <=5) 
	group3 = (5 < age) & (age <10) 
	group4 = age >= 10  
	
	'''
	keep if not smoothing
	star1_r = star_r[group1]
	star1_vrot = star_vrot[group1]
	
	star2_r = star_r[group2]
	star2_vrot = star_vrot[group2]
	
	star3_r = star_r[group3]
	star3_vrot = star_vrot[group3]
	
	star4_r = star_r[group4]
	star4_vrot = star_vrot[group4]
	'''
	
	#get smoothed data
	star1vx, star1vy, star1vz, star1x, star1y, star1z, star1_disp_Z, star1_disp_R, star1_disp_phi = smoothing(group1, 275)
	star2vx, star2vy, star2vz, star2x, star2y, star2z, star2_disp_Z, star2_disp_R, star2_disp_phi = smoothing(group2, 200)
	star3vx, star3vy, star3vz, star3x, star3y, star3z, star3_disp_Z, star3_disp_R, star3_disp_phi = smoothing(group3, 200)
	star4vx, star4vy, star4vz, star4x, star4y, star4z, star4_disp_Z, star4_disp_R, star4_disp_phi = smoothing(group4, 275)

	np.savetxt('/Volumes/Titan/analogs/TNGdata/dispersion/{}_star1_particle_dispersion.txt'.format(halo), np.c_[star1x, star1y, star1_disp_R,star1_disp_Z, star1_disp_phi,], fmt='%1.16f', delimiter=' ', header='x, y, radial component of dispersion, z component LOS), azimuthal component')
	np.savetxt('/Volumes/Titan/analogs/TNGdata/dispersion/{}_star2_particle_dispersion.txt'.format(halo), np.c_[star2x, star2y, star2_disp_R,star2_disp_Z, star2_disp_phi,], fmt='%1.16f', delimiter=' ', header='x, y, radial component of dispersion, z component LOS), azimuthal component')
	np.savetxt('/Volumes/Titan/analogs/TNGdata/dispersion/{}_star3_particle_dispersion.txt'.format(halo), np.c_[star3x, star3y, star3_disp_R,star3_disp_Z, star3_disp_phi,], fmt='%1.16f', delimiter=' ', header='x, y, radial component of dispersion, z component LOS), azimuthal component')
	np.savetxt('/Volumes/Titan/analogs/TNGdata/dispersion/{}_star4_particle_dispersion.txt'.format(halo), np.c_[star4x, star4y, star4_disp_R,star4_disp_Z, star4_disp_phi,], fmt='%1.16f', delimiter=' ', header='x, y, radial component of dispersion, z component LOS), azimuthal component')

	gasvx, gasvy, gasvz, gasx, gasy, gasz, gas_disp = gas_smoothing(275) 

	star1_med_disp_Z.append(np.nanmedian(star1_disp_Z))
	star2_med_disp_Z.append(np.nanmedian(star2_disp_Z))
	star3_med_disp_Z.append(np.nanmedian(star3_disp_Z))
	star4_med_disp_Z.append(np.nanmedian(star4_disp_Z))

	star1_med_disp_R.append(np.nanmedian(star1_disp_R))
	star2_med_disp_R.append(np.nanmedian(star2_disp_R))
	star3_med_disp_R.append(np.nanmedian(star3_disp_R))
	star4_med_disp_R.append(np.nanmedian(star4_disp_R))

	star1_med_disp_phi.append(np.nanmedian(star1_disp_phi))
	star2_med_disp_phi.append(np.nanmedian(star2_disp_phi))
	star3_med_disp_phi.append(np.nanmedian(star3_disp_phi))
	star4_med_disp_phi.append(np.nanmedian(star4_disp_phi))

	#smoothed 
	star1_r = distance(star1x / h, star1y / h, star1z / h)
	star1_vrad_smoothed = radial_v(star1x, star1y, star1z, star1vx, star1vy, star1vz)
	star1_vrot = vrot(star1x, star1y, star1z, star1vx, star1vy, star1vz)
	
	star2_r = distance(star2x / h, star2y / h, star2z / h)
	star2_vrad_smoothed = radial_v(star2x, star2y, star2z, star2vx, star2vy, star2vz)
	star2_vrot = vrot(star2x, star2y, star2z, star2vx, star2vy, star2vz)
	
	star3_r = distance(star3x / h, star3y / h, star3z / h)
	star3_vrad_smoothed = radial_v(star3x, star3y, star3z, star3vx, star3vy, star3vz)
	star3_vrot = vrot(star3x, star3y, star3z, star3vx, star3vy, star3vz)
	
	star4_r = distance(star4x / h, star4y / h, star4z / h)
	star4_vrad_smoothed = radial_v(star4x, star4y, star4z, star4vx, star4vy, star4vz)
	star4_vrot = vrot(star4x, star4y, star4z, star4vx, star4vy, star4vz)
	
	gas_smoothed_r = distance(gasx / h, gasy / h, gasz / h)
	gas_vrad_smoothed = radial_v(gasx, gasy, gasz, gasvx, gasvy, gasvz)
	gas_smoothed_vrot = vrot(gasx, gasy, gasz, gasvx, gasvy, gasvz)

	#divide stars and gas into radial and average the rotation velocities within each bin 
	R_min = 2 #kpc
	R_max = 21
	delta_r= 0.1 #kpc
	
	r_bins=np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)
	
	#will contain the average rotation velocity at each radial bin
	star1_avg_vrot = np.zeros(len(r_bins))
	star2_avg_vrot = np.zeros(len(r_bins))
	star3_avg_vrot = np.zeros(len(r_bins))
	star4_avg_vrot = np.zeros(len(r_bins))
	gas_avg_vrot = np.zeros(len(r_bins))
	gas_avg_vrot_smoothed = np.zeros(len(r_bins))
	
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
		gas_vrots_smoothed=[b for b, a in zip(gas_smoothed_vrot, gas_smoothed_r) if a>=r_bins[i] and a<r_bins[i+1]]
		gas_avg_vrot_smoothed[i]=np.median(gas_vrots_smoothed)
	
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

	#using smoothed gas
	star1_ad_smoothed = va(gas_avg_vrot_smoothed, star1_avg_vrot)
	star2_ad_smoothed = va(gas_avg_vrot_smoothed, star2_avg_vrot)
	star3_ad_smoothed = va(gas_avg_vrot_smoothed, star3_avg_vrot)
	star4_ad_smoothed = va(gas_avg_vrot_smoothed, star4_avg_vrot)
	
	#remove nan values (that resulted from having empty radial bins)
	star1_ad_smoothed = star1_ad_smoothed[~np.isnan(star1_ad_smoothed)]
	star2_ad_smoothed = star2_ad_smoothed[~np.isnan(star2_ad_smoothed)]
	star3_ad_smoothed = star3_ad_smoothed[~np.isnan(star3_ad_smoothed)]
	star4_ad_smoothed = star4_ad_smoothed[~np.isnan(star4_ad_smoothed)] 
	
	#plots-- 1 spatial map and 1 plot of rotation curves and 1 histogram of AD
	single_plot()
	plt.scatter(r_bins[:-1], gas_avg_vrot_smoothed[:-1], c = 'darkgrey', s=12, label='gas')
	plt.scatter(r_bins[:-1], star1_avg_vrot[:-1], c = 'b', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
	plt.scatter(r_bins[:-1], star2_avg_vrot[:-1], c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
	plt.scatter(r_bins[:-1], star3_avg_vrot[:-1], c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
	plt.scatter(r_bins[:-1], star4_avg_vrot[:-1], c = 'r', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
	plt.ylim(0,250)
	plt.xlim(0,22)
	plt.ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=16)
	plt.xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=16)
	plt.legend(loc=2, frameon=False, fontsize=11)
	plt.savefig('/Volumes/Titan/analogs/TNGdata/plots/rcs/{}_rc.png'.format(int(halo)), bbox_inches='tight') #/Volumes/FRIEND/analogs/plots/TNG/rcs/
	plt.close()
	
	#========================================================
	single_plot()
	plt.hist(star1_ad_smoothed[:-1], bins=range(-200, 300, 20), label='Group 1: {} stars, AD ={}'.format(sum(group1), np.median(star1_ad_smoothed[:-1])), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
	plt.hist(star2_ad_smoothed[:-1], bins=range(-200, 300, 20), label='Group 2: {} stars, AD ={}'.format(sum(group2), np.median(star2_ad_smoothed[:-1])), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
	plt.hist(star3_ad_smoothed[:-1], bins=range(-200, 300, 20), label='Group 3: {} stars, AD ={}'.format(sum(group3), np.median(star3_ad_smoothed[:-1])), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
	plt.hist(star4_ad_smoothed[:-1], bins=range(-200, 300, 20), label='Group 4: {} stars, AD ={}'.format(sum(group4), np.median(star4_ad_smoothed[:-1])), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
	plt.legend(loc=2, frameon=False)
	plt.xlim(-300,300)
	plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
	plt.savefig('/Volumes/Titan/analogs/TNGdata/plots/ad/{}_AD.png'.format(int(halo)), bbox_inches='tight')

	#combinging rc and ad as one plot
	f, axes= plt.subplots(1,2, sharey=False, sharex=False, figsize=(18,6))
	for ax in axes:
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
	    	    ax.spines[axis].set_linewidth(1)
	axes[0].scatter(r_bins[:-1], gas_avg_vrot_smoothed[:-1], c = 'darkgrey', s=12, label='gas')
	axes[0].scatter(r_bins[:-1], star1_avg_vrot[:-1], c = 'b', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
	axes[0].scatter(r_bins[:-1], star2_avg_vrot[:-1], c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
	axes[0].scatter(r_bins[:-1], star3_avg_vrot[:-1], c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
	axes[0].scatter(r_bins[:-1], star4_avg_vrot[:-1], c = 'r', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
	axes[0].set_ylim(0,200)
	axes[0].set_xlim(0,20)
	axes[0].set_ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=16)
	axes[0].set_xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=16)
	axes[0].legend(loc=2, frameon=False, fontsize=13)

	axes[1].hist(star1_ad_smoothed[:-1], bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_ad_smoothed[:-1]),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
	axes[1].hist(star2_ad_smoothed[:-1], bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_ad_smoothed[:-1]),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
	axes[1].hist(star3_ad_smoothed[:-1], bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_ad_smoothed[:-1]),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
	axes[1].hist(star4_ad_smoothed[:-1], bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_ad_smoothed[:-1]),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
	axes[1].legend(loc=1, frameon=False, fontsize=1)
	axes[1].set_xlim(-115, 150)
	axes[1].set_xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=16)

	#plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('/Volumes/Titan/analogs/TNGdata/plots/rc_ad/{}_AD_RC.png'.format(int(halo)), bbox_inches='tight')
	plt.close()

	#======================================================
	rc('font', family = 'serif')
	f, axes= plt.subplots(3,4, sharey=True, sharex=False, figsize=(14,10.4))
	
	centerx= 17
	centery= -16
	radius1= 200 / 60 / 60 * 13.86 
	radius2= 275 / 60 / 60 * 13.86 
	c = plt.Circle((centerx, centery), radius2, color='k', fill=False)
	c1 = plt.Circle((centerx, centery), radius1, color='k', fill=False)
	c2 = plt.Circle((centerx, centery), radius2, color='k', fill=False)
	c3 = plt.Circle((centerx, centery), radius1, color='k', fill=False)
	c4 = plt.Circle((centerx, centery), radius1, color='k', fill=False)
	c5 = plt.Circle((centerx, centery), radius1, color='k', fill=False)
	c6 = plt.Circle((centerx, centery), radius2, color='k', fill=False)
	c7 = plt.Circle((centerx, centery), radius2, color='k', fill=False)
	
	axes[0,0].scatter(star_x[group1] / h, star_y[group1] / h, c=star_vz[group1], cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[0,1].scatter(star_x[group2] / h, star_y[group2] / h, c=star_vz[group2], cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[0,2].scatter(star_x[group3] / h, star_y[group3] / h, c=star_vz[group3], cmap='plasma', s=4, vmin=-150,vmax=100) 
	im0=axes[0,3].scatter(star_x[group4] / h, star_y[group4] / h, c=star_vz[group4], cmap='plasma', s=4, vmin=-150,vmax=100)
	
	#axes[0,0].annotate('Group 1', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[0,1].annotate('Group 2', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[0,2].annotate('Group 3', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[0,3].annotate('Group 4', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	
	axes[1,0].add_artist(c)
	axes[1,0].scatter(star1x / h, star1y / h, c= star1vz, cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[1,1].scatter(star2x / h, star2y / h, c= star2vz, cmap='plasma', s=4, vmin=-150,vmax=100)
	axes[1,1].add_artist(c4)
	axes[1,2].scatter(star3x / h, star3y / h, c= star3vz, cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[1,2].add_artist(c5)
	im1=axes[1,3].scatter(star4x / h, star4y / h, c= star4vz, cmap='plasma', s=4, vmin=-150,vmax=100)
	axes[1,3].add_artist(c6)
	
	#axes[1,0].annotate('Group 1', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[1,1].annotate('Group 2', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[1,2].annotate('Group 3', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[1,3].annotate('Group 4', xy=(-19.5,16), horizontalalignment='left', fontsize=12)

	#to look at weird LOS rotation pattern
	axes[2,0].scatter(star1x / h, star1y / h, c= star1_vrot, cmap='viridis', s=4, vmin=-0,vmax=300) 
	axes[2,1].scatter(star2x / h, star2y / h, c= star2_vrot, cmap='viridis', s=4, vmin=-0,vmax=300) 
	axes[2,2].scatter(star3x / h, star3y / h, c= star3_vrot, cmap='viridis', s=4, vmin=-0,vmax=300) 
	im2=axes[2,3].scatter(star4x / h, star4y / h, c= star4_vrot, cmap='viridis', s=4, vmin=-0,vmax=300)
	
	#axes[2,0].annotate('Group 1',xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[2,1].annotate('Group 2', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[2,2].annotate('Group 3', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	#axes[2,3].annotate('Group 4', xy=(-19.5,16), horizontalalignment='left', fontsize=12)

	axes[0,0].set_title(r'$\rm Group\ 1: \leq 1\ Gyr$')
	axes[0,1].set_title(r'$\rm Group\ 2: 1-5\ Gyr$')
	axes[0,2].set_title(r'$\rm Group\ 3: 5-10\ Gyr $')
	axes[0,3].set_title(r'$\rm Group\ 4: \geq 10\ Gyr $')
	
	for ax in axes[0,:]:
		ax.set_xlim(20, -20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
		ax.set_xticklabels([])
	
	for ax in axes[1,:]:
		ax.set_xlim(20, -20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
		ax.set_xticklabels([])
	
	for ax in axes[2,:]:
		ax.set_xlim(20, -20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		nbins = 5 
		ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	for ax in axes[:,0]:
		ax.set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
		ax.set_ylim(-20,20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		nbins = 7
		ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	for ax in axes[:,1]:
		ax.set_ylim(-20,20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	for ax in axes[:,2]:
		ax.set_ylim(-20,20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	for ax in axes[:,3]:
		ax.set_ylim(-20,20)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
		ax.tick_params(axis='y',which='both',right='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	f.subplots_adjust(right=0.885)
	cbar_ax1 = f.add_axes([0.89,0.38,0.015,0.512])
	cbar_ax2 = f.add_axes([0.89,0.115,0.015,0.25])
	clb1=f.colorbar(im1, cax=cbar_ax1)
	clb2=f.colorbar(im2, cax=cbar_ax2)
	clb1.set_label(r'$\rm Individual,\ Median\ LOS\ velocity:\ v, \ \overline{v}\ (km\ s^{-1})$', fontsize=13)
	clb2.set_label(r'$\rm Rotation\ Velocity: \ v_{rot}\ (km\ s^{-1})$', fontsize=13, labelpad=12)
	axes[0,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[0,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[0,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[0,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2,0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2,1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2,2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2,3].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('/Volumes/Titan/analogs/TNGdata/plots/maps/{}_star_map.png'.format(int(halo)), bbox_inches='tight')

	#========================================================================================
	#examine the gas
	rc('font', family = 'serif')
	f, axes= plt.subplots(1,3, sharey=True, sharex=False)#, figsize=(14,10.4))

	im0 = axes[0].scatter(gas_x / h, gas_y / h, c=gas_vz, cmap='plasma', s=4, vmin=-150,vmax=100)
	axes[1].scatter(gasx / h, gasy / h, c=gasvz, cmap='plasma', s=4, vmin=-150,vmax=100)
	im1 = axes[2].scatter(gasx / h, gasy / h, c=gas_smoothed_vrot, cmap='viridis', s=4, vmin=-0,vmax=330) 

	axes[0].annotate('Individual',xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	axes[1].annotate('Smoothed', xy=(-19.5,16), horizontalalignment='left', fontsize=12)
	axes[2].annotate('Rotation', xy=(-19.5,16), horizontalalignment='left', fontsize=12)

	axes[0].set_ylim(-18, 18)
	axes[0].set_xlim(18, -18)
	axes[0].set(adjustable='box-forced', aspect='equal')
	axes[0].tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	axes[0].tick_params(axis='x',which='both',top='on', direction='in')
	axes[0].tick_params(which='both', width=1)
	axes[0].tick_params(which='major', length=7)
	axes[0].tick_params(which='minor', length=4)
	axes[0].tick_params(labelsize=12) 
	axes[0].minorticks_on()
	
	for ax in axes[:]:
		ax.set_xlim(18, -18)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)

	f.subplots_adjust(right=0.885)
	cbar_ax1 = f.add_axes([0.89,0.375,0.015,0.512])
	cbar_ax2 = f.add_axes([0.89,0.1,0.015,0.27])
	clb1=f.colorbar(im0, cax=cbar_ax1)
	clb2=f.colorbar(im1, cax=cbar_ax2)
	clb1.set_label(r'$\rm LOS\ velocity:\ v, \ \overline{v}\ (km\ s^{-1})$', fontsize=13)
	clb2.set_label(r'$\rm Rotation\ Velocity: \ v_{rot}\ (km\ s^{-1})$', fontsize=13, labelpad=12)
	axes[0].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[1].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	axes[2].scatter(0,0,marker='+', c='dodgerblue', linewidth=2)
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('/Volumes/Titan/analogs/TNGdata/plots/maps/{}_gas_map.png'.format(int(halo)), bbox_inches='tight')

	# #========================================================================================
	# #to examine the gas on nan ad analogs
	# single_plot()
	# plt.scatter(gas_x / h, gas_y / h, c=gas_vrot, alpha = .5, s=6, vmin=100, vmax=350)
	# clb=plt.colorbar()
	# clb.set_label('km/s')
	# plt.xlim(-20,20)
	# plt.ylim(-20, 20)
	# plt.xlabel('kpc')
	# plt.ylabel('kpc')
	# plt.savefig('/Volumes/FRIEND/analogs/plots/disk_splitting/{}_gas_map.png'.format(int(halo)), bbox_inches='tight')
	# plt.close()
	
	# single_plot()
	# plt.hist(gas_fraction)
	# plt.savefig('/Volumes/FRIEND/analogs/plots/disk_splitting/{}_gas_hist.png'.format(int(halo)), bbox_inches='tight')
	# plt.close()
	# #======================================================================================

	#save data to file
	np.savetxt('/Volumes/Titan/analogs/TNGdata/smoothed_vrot/{}_vrot_smoothed.txt'.format(int(halo)), np.c_[r_bins[:-1], gas_avg_vrot_smoothed[:-1], star1_avg_vrot[:-1], star2_avg_vrot[:-1], star3_avg_vrot[:-1], star4_avg_vrot[:-1]], fmt='%1.16f', delimiter=' ', header='r bin (kpc), avg v_rot gas (km/s), avg v_rot star 1 (km/s), avg v_rot star 2 (km/s), avg v_rot star 3 (km/s), avg v_rot star 4 (km/s)')

np.savetxt('/Volumes/Titan/analogs/TNGdata/dispersion/star_particle_dispersion.txt', np.c_[halos, star1_med_disp_R, star2_med_disp_R, star3_med_disp_R, star4_med_disp_R, star1_med_disp_Z, star2_med_disp_Z, star3_med_disp_Z, star4_med_disp_Z, star1_med_disp_phi, star2_med_disp_phi, star3_med_disp_phi, star4_med_disp_phi], fmt='%1.16f', delimiter=' ', header='halo, radial component of dispersion for stellar group 1 2 3 4, z component, azimuthal component')


