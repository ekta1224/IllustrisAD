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

halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt')


for halo in halos:
	print(halo)
	#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
	star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass_all, star_factor = np.loadtxt('../data/M31analog_{}_star_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True) 
	gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass_all, gas_fraction = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(int(halo)), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True)
	
	#coordinates ckpc/h -- since at z=0, kpc/h
	#speeds km sqrt(a)/s -- since at z=0, km/s
	
	#Illustric value
	h = 0.704
	
	#calculate distance from center-- we don't need to deproject the position (yay simulations)
	def distance(x, y, z):
		#(coordinates are already shifted to 0,0,0)
		return np.sqrt((x**2) + (y**2) + (z**2)) #whatever units the input is in
	
	#calculate the radial velocity
	def radial_v(x, y, z, vx, vy, vz):
		#convert kpc/h to km
		x_km = x * 3.086e+16 / h
		y_km = y * 3.086e+16 / h
		z_km = z * 3.086e+16 / h
		dot_product = (x_km * vx) + (y_km * vy)+ (z_km * vz) 
		r_mag = distance(x_km, y_km, z_km)
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
	close = star_z / h <= 10#) & (star_y / h < 1 - 17 * (star_x / h) / 14) #second part is used for disk splitting tests
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
	neutral_gas = gas_fraction > 0.95 #what is the best limit?
	
	gas_x = gas_x[neutral_gas]
	gas_y = gas_y[neutral_gas]
	gas_z = gas_z[neutral_gas]
	gas_vx = gas_vx[neutral_gas]
	gas_vy = gas_vy[neutral_gas]
	gas_vz = gas_vz[neutral_gas]

	#get rid of particles that are too distant in the z direction
	close = gas_z / h <= 10 #) & (gas_y / h < 1 - 17 * (gas_x / h) / 14) #second part is used for disk splitting tests
	gas_x = gas_x[close]
	gas_y = gas_y[close]
	gas_z = gas_z[close]
	gas_vx = gas_vx[close]
	gas_vy = gas_vy[close]
	gas_vz = gas_vz[close]
	
	gas_vrad = radial_v(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz) #km/s
	gas_r = distance(gas_x / h, gas_y / h, gas_z / h) #kpc
	
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
		star_dispersion = []
		c = SkyCoord(ra = star_xs / 13.86 / h, dec = star_ys / 13.86 / h, unit=(u.deg,u.deg))
		for i in range(len(star_xs)):
			c1 = SkyCoord(star_xs[i] / 13.86 / h, star_ys[i] / 13.86 / h, unit=(u.deg,u.deg)) 
			sep = c1.separation(c)
			good = sep.arcsecond <= circle_size #put stars into smoothing circle of this size
			star_x_velocities = star_vxs[good]
			star_y_velocities = star_vys[good]
			star_z_velocities = star_vzs[good]
			radial_vels = radial_v(star_xs[good], star_ys[good], star_zs[good], star_x_velocities, star_y_velocities, star_z_velocities) #to get dispersion
			if len(star_x_velocities) >= 15:
				star_goodcenter_x.append(star_xs[i]) #kpc/h
				star_goodcenter_y.append(star_ys[i]) #kpc/h
				star_goodcenter_z.append(star_zs[i]) #kpc/h
				star_smoothed_vx.append(np.median(star_x_velocities)) #km/s
				star_smoothed_vy.append(np.median(star_y_velocities)) #km/s
				star_smoothed_vz.append(np.median(star_z_velocities)) #km/s
				star_dispersion.append(np.std(radial_vels))
	
		return np.array((star_smoothed_vx)), np.array((star_smoothed_vy)), np.array((star_smoothed_vz)), np.array((star_goodcenter_x)), np.array((star_goodcenter_y)), np.array((star_goodcenter_z)), np.array((star_dispersion))
	
	#=========================================================================
	
	#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
	def vrot(x, y, z, vx, vy, vz):
		v_rad = radial_v(x, y, z, vx, vy, vz) #km/s
		v_tot = np.sqrt((vx**2) + (vy**2) + (vz**2)) #km/s
		v_rot = np.sqrt((v_tot**2) - (v_rad**2))
		return v_rot
	
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
	
	#divide the stars into 4 age bins
	group1 = age <= 1 #Gyr
	group2 = (1 <=age) & (age <=5) 
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
	star1vx, star1vy, star1vz, star1x, star1y, star1z, star1_disp = smoothing(group1, 275)
	star2vx, star2vy, star2vz, star2x, star2y, star2z, star2_disp = smoothing(group2, 200)
	star3vx, star3vy, star3vz, star3x, star3y, star3z, star3_disp = smoothing(group3, 200)
	star4vx, star4vy, star4vz, star4x, star4y, star4z, star4_disp = smoothing(group4, 275)
	
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
	
	#divide stars and gas into radial and average the rotation velocities within each bin 
	R_min = 0 #kpc
	R_max = 21
	delta_r= 0.1 #kpc
	
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
	
	#plots-- 1 spatial map and 1 plot of rotation curves and 1 histogram of AD
	single_plot()
	plt.scatter(r_bins[:-1], gas_avg_vrot[:-1], c = 'darkgrey', s=12, label='gas')
	plt.scatter(r_bins[:-1], star1_avg_vrot[:-1], c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1')
	plt.scatter(r_bins[:-1], star2_avg_vrot[:-1], c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2')
	plt.scatter(r_bins[:-1], star3_avg_vrot[:-1], c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3')
	plt.scatter(r_bins[:-1], star4_avg_vrot[:-1], c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4')
	plt.ylim(0,500)
	plt.xlim(0,20)
	plt.ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=13)
	plt.xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
	plt.legend(loc=2, frameon=False)
	plt.savefig('/Volumes/FRIEND/analogs/plots/rcs/{}_rc.png'.format(int(halo)), bbox_inches='tight')
	plt.close()
	
	#========================================================
	single_plot()
	plt.hist(star1_ad[:-1], bins=range(-200, 300, 20), label='Group 1: {} stars, AD ={}'.format(sum(group1), np.median(star1_ad[:-1])), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
	plt.hist(star2_ad[:-1], bins=range(-200, 300, 20), label='Group 2: {} stars, AD ={}'.format(sum(group2), np.median(star2_ad[:-1])), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
	plt.hist(star3_ad[:-1], bins=range(-200, 300, 20), label='Group 3: {} stars, AD ={}'.format(sum(group3), np.median(star3_ad[:-1])), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
	plt.hist(star4_ad[:-1], bins=range(-200, 300, 20), label='Group 4: {} stars, AD ={}'.format(sum(group4), np.median(star4_ad[:-1])), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
	plt.legend(loc=2, frameon=False)
	plt.xlim(-300,300)
	plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
	plt.savefig('/Volumes/FRIEND/analogs/plots/hists/{}_AD.png'.format(int(halo)), bbox_inches='tight')
	
	#======================================================
	rc('font', family = 'serif')
	f, axes= plt.subplots(3,4, sharey=True, sharex=False, figsize=(14,10.4))
	
	centerx= 13
	centery= -12
	radius1= 200 / 60 / 60 * 13.86 #fix scaling factor
	radius2= 275 / 60 / 60 * 13.86 #fix scaling factor
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
	
	axes[0,0].annotate('Group 1', xy=(-.55,13.5), horizontalalignment='right', fontsize=12)
	axes[0,1].annotate('Group 2', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[0,2].annotate('Group 3', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[0,3].annotate('Group 4', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	
	axes[1,0].add_artist(c)
	axes[1,0].scatter(star1x / h, star1y / h, c= star1vz, cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[1,1].scatter(star2x / h, star2y / h, c= star2vz, cmap='plasma', s=4, vmin=-150,vmax=100)
	axes[1,1].add_artist(c4)
	axes[1,2].scatter(star3x / h, star3y / h, c= star3vz, cmap='plasma', s=4, vmin=-150,vmax=100) 
	axes[1,2].add_artist(c5)
	im1=axes[1,3].scatter(star4x / h, star4y / h, c= star4vz, cmap='plasma', s=4, vmin=-150,vmax=100)
	axes[1,3].add_artist(c6)
	
	axes[1,0].annotate('Group 1', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[1,1].annotate('Group 2', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[1,2].annotate('Group 3', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[1,3].annotate('Group 4', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	
	axes[2,0].scatter(star1x / h, star1y / h, c= star1_disp, cmap='copper', s=4, vmin=0,vmax=200) 
	axes[2,0].add_artist(c2)
	axes[2,1].scatter(star2x / h, star2y / h, c= star2_disp, cmap='copper', s=4, vmin=0,vmax=200)
	axes[2,1].add_artist(c1)
	axes[2,2].add_artist(c3)
	axes[2,2].scatter(star3x / h, star3y / h, c= star3_disp, cmap='copper', s=4,vmin=0,vmax=200) 
	axes[2,3].add_artist(c7)
	im2=axes[2,3].scatter(star4x / h, star4y / h, c= star4_disp, cmap='copper', s=4,vmin=0,vmax=200)
	
	axes[2,0].annotate('Group 1',xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[2,1].annotate('Group 2', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[2,2].annotate('Group 3', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	axes[2,3].annotate('Group 4', xy=(-.550,13.5), horizontalalignment='right', fontsize=12)
	
	for ax in axes[0,:]:
		ax.set_xlim(18, -18)
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
	
	for ax in axes[1,:]:
		ax.set_xlim(18, -18)
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
	
	for ax in axes[2,:]:
		ax.set_xlim(18, -18)
		ax.set(adjustable='box-forced', aspect='equal')
		ax.set_xlabel(r'$\xi\ (kpc)$', fontsize=13)
		ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
		ax.tick_params(axis='x',which='both',top='on', direction='in')
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)
		ax.tick_params(labelsize=12) 
		ax.minorticks_on()
		for axis in ['top','bottom','left','right']:
		        ax.spines[axis].set_linewidth(1)
	
	for ax in axes[:,0]:
		ax.set_ylabel(r'$\eta\ (kpc)$', fontsize=13)
		ax.set_ylim(-18,18)
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
		ax.set_ylim(-18,18)
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
		ax.set_ylim(-18,18)
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
		ax.set_ylim(-18,18)
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
	cbar_ax1 = f.add_axes([0.89,0.375,0.015,0.512])
	cbar_ax2 = f.add_axes([0.89,0.1,0.015,0.27])
	clb1=f.colorbar(im1, cax=cbar_ax1)
	clb2=f.colorbar(im2, cax=cbar_ax2)
	clb1.set_label(r'$\rm Individual,\ Mean\ LOS\ velocity:\ v, \ \overline{v}\ (km\ s^{-1})$', fontsize=13)
	clb2.set_label(r'$\rm Velocity\ Dispersion: \sigma\ (km\ s^{-1})$', fontsize=13, labelpad=12)
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
	plt.savefig('/Volumes/FRIEND/analogs/plots/maps/{}_map.png'.format(int(halo)), bbox_inches='tight')

	#========================================================================================
	#to examine the gas on nan ad analogs
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
	#======================================================================================

	#save data to file
	np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot.txt'.format(int(halo)), np.c_[r_bins[:-1], gas_avg_vrot[:-1], star1_avg_vrot[:-1], star2_avg_vrot[:-1], star3_avg_vrot[:-1], star4_avg_vrot[:-1]], fmt='%1.16f', delimiter=' ', header='r bin (kpc), avg v_rot gas (km/s), avg v_rot star 1 (km/s), avg v_rot star 2 (km/s), avg v_rot star 3 (km/s), avg v_rot star 4 (km/s)')


