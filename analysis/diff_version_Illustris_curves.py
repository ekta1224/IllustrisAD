import numpy as np 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 
from astropy.coordinates import SkyCoord
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

'''
-divides stellar particles into 4 age bins
-eliminates gas particles that don't have a 'high' fraction of neutral hydrogen
-construct a rotation curve for the gas and star particles, binned radially
-calculates the asymmetric drift at each radial pin and outputs a histogram to compare the four age bins
'''

# Hello amanda take a deep breath, debugging starts now!! 

halo = '361428'

#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass_all, star_factor = np.loadtxt('../data/M31analog_{}_star_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True) 
gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass_all, gas_fraction = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True)

#coordinates ckpc/h -- since at z=0, kpc/h
#speeds km sqrt(a)/s -- since at z=0, km/s

#Illustric value
h = 0.704

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y):#, z):
	#convert to kpc (coordinates are already shifted to 0,0,0)
	#z_kpc = z / h
	return np.sqrt((x**2) + (y**2)) #+ (z_kpc**2)) #whatever units x and y are in 

#calculate the radial velocity
def radial_v(x, y, vx, vy):
	#convert kpc/h to km
	x_km = x * 3.086e+16 / h
	y_km = y * 3.086e+16 / h
	dot_product = (x_km * vx) + (y_km * vy)# + (z * vz) 
	r_mag = distance(x_km, y_km)#, z)
	return dot_product / r_mag #km/s

# ENIA: particles are willddddddd 

#get rid of star partiles that are really wind
not_wind = star_factor > 0
star_x = star_x[not_wind]
star_y = star_y[not_wind]
#star_z = star_z[not_wind]
star_vx = star_vx[not_wind]
star_vy = star_vy[not_wind]
#star_vz = star_vz[not_wind]
star_factor = star_factor[not_wind]

#star_vrad = radial_v(star_x, star_y, star_z, star_vx, star_vy, star_vz)

#star_r = distance(star_x, star_y, star_z)

#only want particles/cells that correspond to a high HI gas fraction
neutral_gas = gas_fraction > 0.95 #what is the best limit?

# ENIA: You are beautiful in every single way words can't even explain 

gas_x = gas_x[neutral_gas]
gas_y = gas_y[neutral_gas]
#gas_z = gas_z[neutral_gas]
#gas_r = distance(gas_x, gas_y, gas_z)
gas_vx = gas_vx[neutral_gas]
gas_vy = gas_vy[neutral_gas]
#gas_vz = gas_vz[neutral_gas]
#gas_vrad = radial_v(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz)

#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
def vrot(x, y, vx, vy):
	v_rad = radial_v(x, y, vx, vy) #km/s
	v_tot = np.sqrt((vx**2) + (vy**2))# + (vz**2))
	v_rot = np.sqrt((v_tot**2) - (v_rad**2))
	return v_rot #km/s

def vrot_with_smoothed_vrad(smoothed_vrad, smoothed_vtot): #used to imitate getting rotation velocity from smoothed radial velocity
	v_rot = np.sqrt((smoothed_vtot**2) - (smoothed_vrad**2))
	return v_rot

# star_vrot = vrot(star_x, star_y, star_z, star_vx, star_vy, star_vz)
# gas_vrot = vrot(gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz)

# ENIA: I belive you can fly 

#calculate the ages of stars from the scale factor
def star_age(scale_factor):
	zform = [(1./a) -1. for a in scale_factor]
	cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
	tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
	tage = np.array([13.8-t for t in tform])
	return tage

age = star_age(star_factor)

#divide the stars into 4 age bins -- do we want larger bins? -- yes do this later amanda
group1 = (0 <= age) & (age <=.06) #average age 30 Myr 
group2 = (.1 <= age) & (age <=.9) #average age 400 Myr
group3 = (1 <= age) & (age <= 3) #average age 2 Gyr
group4 = (3 <= age) & (age <= 11) #average age 4 Gyr

# star1_r = star_r[group1]
# star1_vrot = star_vrot[group1]

# star2_r = star_r[group2]
# star2_vrot = star_vrot[group2]

# star3_r = star_r[group3]
# star3_vrot = star_vrot[group3]

# star4_r = star_r[group4]
# star4_vrot = star_vrot[group4]

#star1_x = star_x[group1]
#star1_y = star_y[group1]
#star1_z = star_z[group1]
#star1_vx = star_vx[group1]
#star1_vy = star_vy[group1]
#star1_vz = star_vz[group1]

# ENIA: the most magical things start with a single click 

#pair gas and star particles 
def pair_particles(starx, stary, gasx, gasy):
	#position ckpc/h -> deg; this isn't exactly right; value is for M31
	star_coords = SkyCoord(ra = starx / 13.67 / h, dec = stary / 13.67 / h, unit=(u.deg,u.deg))
	gas_coords = SkyCoord(ra = gasx / 13.67 / h, dec = gasy / 13.67 / h, unit=(u.deg,u.deg))
	idx, d2d, d3d = star_coords.match_to_catalog_sky(gas_coords)
	#make sure paired points aren't too far away
	return idx, d2d

idx1, d2d1 = pair_particles(star_x[group1], star_y[group1], gas_x, gas_y)
idx2, d2d2 = pair_particles(star_x[group2], star_y[group2], gas_x, gas_y)
idx3, d2d3 = pair_particles(star_x[group3], star_y[group3], gas_x, gas_y)
idx4, d2d4 = pair_particles(star_x[group4], star_y[group4], gas_x, gas_y)

def smoothing(group, index, pair_distance): #we do not smooth the HI
	star_xs = star_x[group]
	star_ys = star_y[group]
	star_vxs = star_vx[group]
	star_vys = star_vy[group]
	gas_vxs = gas_vx[index]
	gas_vys = gas_vy[index]

	#only want to look at pairs that are actually close to each other
	close_pairs = pair_distance < 40 #actually think of a good number to do; seems like going from 50 -> 20 does not help the shape

	star_xs = star_xs[close_pairs] 
	star_ys = star_ys[close_pairs] 
	star_vxs =star_vxs[close_pairs]
	star_vys =star_vys[close_pairs]
	gas_vxs = gas_vxs[close_pairs] 
	gas_vys = gas_vys[close_pairs] 

	star_smoothed_vrad = []
	star_smoothed_vrot = []
	gas_vrad = []
	gas_vrot = []
	star_x_goodcenter = []
	star_y_goodcenter = []
	star_dispersion = []
	star_r = []
	c = SkyCoord(ra = star_xs / 13.67 / h, dec = star_ys / 13.67 / h, unit=(u.deg,u.deg))
	for i in range(len(star_xs)):
		c1 = SkyCoord(star_xs[i] / 13.67 / h, star_ys[i] / 13.67 / h, unit=(u.deg,u.deg)) #go through all coordinates one at a time
		sep = c1.separation(c)
		good = sep.arcsecond < 200 #put stars into smoothing circle of this size
		star_x_velocities = star_vxs[good]
		star_y_velocities = star_vys[good]
		star_radial_vels = radial_v(star_xs[good], star_ys[good], star_x_velocities, star_y_velocities)
		#star_vtots = np.sqrt((star_x_velocities**2) + (star_y_velocities**2))
		if len(star_x_velocities) > 15: #only want circles with at least 15 stars
			star_x_goodcenter.append(star_xs[i]) #kpc
			star_y_goodcenter.append(star_ys[i]) #kpc
			star_r.append(distance(star_xs[i] / h, star_ys[i] / h))#, star_z[i]))
			star_smoothed_vrad.append(np.median(star_radial_vels)) #average the radial velocity
			#star_smoothed_vrot.append(vrot_with_smoothed_vrad(np.median(star_radial_vels), np.median(star_vtots))) #not smoothing the rotation velocity but using averaged radial and total velocites
			star_smoothed_vrot.append(vrot(star_xs[i], star_ys[i], np.median(star_x_velocities), np.median(star_y_velocities))) #which one is better? these look the same
			star_dispersion.append(np.std(star_radial_vels))
			gas_vrad.append(radial_v(star_xs[i], star_ys[i], gas_vxs[i], gas_vys[i])) #don't smooth gas but keep point with good center
			gas_vrot.append(vrot(star_xs[i], star_ys[i], gas_vxs[i], gas_vys[i])) #don't smooth gas but keep point with good center
		
	return star_x_goodcenter, star_y_goodcenter, star_r, star_smoothed_vrad, star_smoothed_vrot, gas_vrad, gas_vrot,  star_dispersion

group1_data = smoothing(group1, idx1, d2d1.arcsecond)
group2_data = smoothing(group2, idx2, d2d2.arcsecond)
group3_data = smoothing(group3, idx3, d2d3.arcsecond)
group4_data = smoothing(group4, idx4, d2d4.arcsecond)

#make rotation curves
delta_r = 0.5 #kpc
bins = np.arange(5, 20, delta_r)
median_r = bins - delta_r / 2

def median_line(r, v):
	medians = np.zeros_like(bins)
	for i in range(len(bins)):
		data = [b for a,b in zip(r, v) if a < bins[i] and a > bins[i] - delta_r]
		medians[i] = np.median(data)
	return medians 

group1_star_vrot_med = median_line(group1_data[2], group1_data[4])
group2_star_vrot_med = median_line(group2_data[2], group2_data[4])
group3_star_vrot_med = median_line(group3_data[2], group3_data[4])
group4_star_vrot_med = median_line(group4_data[2], group4_data[4])
group1_gas_vrot_med = median_line(group1_data[2], group1_data[6])
group2_gas_vrot_med = median_line(group2_data[2], group2_data[6])
group3_gas_vrot_med = median_line(group3_data[2], group3_data[6])
group4_gas_vrot_med = median_line(group4_data[2], group4_data[6])

rc('font', family = 'serif')
f, axes= plt.subplots(4,1, sharey=False, sharex=True, figsize=(4,9.8))
axes[0].scatter(group1_data[2], group1_data[6], s=2, c='darkgray')
axes[0].scatter(group1_data[2], group1_data[4], s=2, c='b', alpha=0.4)
axes[0].plot(median_r, group1_star_vrot_med, linestyle='-', c='black', linewidth = 1.8, alpha=.85)
axes[0].plot(median_r, group1_gas_vrot_med, linestyle='--', c='black', linewidth = 1.8, alpha=.85)
axes[1].scatter(group2_data[2], group2_data[6], s=2, c='darkgray')
axes[1].scatter(group2_data[2], group2_data[4], s=2, c='m', alpha=0.4)
axes[1].plot(median_r, group2_star_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[1].plot(median_r, group2_gas_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[2].scatter(group3_data[2], group3_data[6], s=2, c='darkgray')
axes[2].scatter(group3_data[2], group3_data[4], s=2, c='green', alpha=0.4)
axes[2].plot(median_r, group3_star_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[2].plot(median_r, group3_gas_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[3].scatter(group4_data[2], group4_data[6], s=2, c='darkgray')
axes[3].scatter(group4_data[2], group4_data[4], s=2, c='r', alpha=0.4)
axes[3].plot(median_r, group4_star_vrot_med, linestyle='-', c='black', linewidth = 1.8)
axes[3].plot(median_r, group4_gas_vrot_med, linestyle='--', c='black', linewidth = 1.8)
axes[0].annotate('30 Myr', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[1].annotate('400 Myr', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[2].annotate('2 Gyr', xy=(19,115), horizontalalignment='right', fontsize=12)
axes[3].annotate('4 Gyr', xy=(19,115), horizontalalignment='right', fontsize=12)


# ENIA: plot, phot, code, trot, mode 

for ax in axes:
	ax.set_xlim(4, 20)
	ax.set_ylim(100,300)
	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	ax.tick_params(which='both', width=2)
	ax.tick_params(which='major', length=7)
	ax.tick_params(which='minor', length=4)
	ax.tick_params(labelsize=12) 
	ax.minorticks_on()
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
axes[3].set_xlabel(r'$\rm Radial\ Distance:\ \it r \ \rm(kpc)$', fontsize=13)
nbins = len(axes[0].get_yticklabels())-1
axes[3].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[1].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
axes[2].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
f.subplots_adjust(left=0.17)
f.text(0.008, 0.5, r'$\rm Rotation\ Velocity:\ \it v_{\rm rot}\ \rm(km\ s^{-1})$', va='center', rotation='vertical', fontsize=13)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/{}_rotation_curve.pdf'.format(halo))
plt.close()

#calculate the asymmetric drift
def va(v_gas, v_star):
	return np.array((v_gas)) - np.array((v_star))

def errors(AD):
	result = np.percentile(AD, [16, 50, 84])
	median = result[1]
	lower_error = result[1] - result[0]
	upper_error = result[2] - result[1]
	return median, lower_error, upper_error

#group1_ad = va(group1_data[4], group1_data[6])
group2_ad = va(group2_data[4], group2_data[6])
group3_ad = va(group3_data[4], group3_data[6])
group4_ad = va(group4_data[4], group4_data[6])

#group1_ad_error = errors(group1_ad)
group2_ad_error = errors(group2_ad)
group3_ad_error = errors(group3_ad)
group4_ad_error = errors(group4_ad)

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
#plt.hist(group1_ad, bins=range(-70, 70, 10), label='30 Myr={} stars'.format(len(group1_ad)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(group2_ad, bins=range(-70, 70, 10), label='400 Myr={} stars, {}'.format(len(group2_ad), round(np.median(group2_ad), 2)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m')
plt.hist(group3_ad, bins=range(-70, 70, 10), label='2 Gyr={} stars, {}'.format(len(group3_ad), round(np.median(group3_ad), 2)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='k')
plt.hist(group4_ad, bins=range(-70, 70, 10), label='4 Gyr={} stars, {}'.format(len(group4_ad), round(np.median(group4_ad), 2)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='r')
plt.legend(loc=1, frameon=False)
plt.xlim(-75,70)
plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/{}_asymmetric_drift_hist.pdf'.format(halo))
plt.close()

#save data to file
#np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot_group1.txt'.format(halo), np.c_[group1_data[0], group1_data[1], group1_data[2], group1_data[3], group1_data[4], group1_data[5], group1_data[6], group1_data[7], #group1_data[8]], fmt='%1.16f', delimiter=' ', header='star_x_goodcenter kpc, star_y_goodcenter, star_r, star_smoothed_vrad km/s, star_smoothed_vrot, gas_vrad, gas_vrot,  star_dispersion')
#np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot_group2.txt'.format(halo), np.c_[group2_data[0], group2_data[1], group2_data[2], group2_data[3], group2_data[4], group2_data[5], group2_data[6], group2_data[7], #group2_data[8]], fmt='%1.16f', delimiter=' ', header='star_x_goodcenter kpc, star_y_goodcenter, star_r, star_smoothed_vrad km/s, star_smoothed_vrot, gas_vrad, gas_vrot,  star_dispersion')
#np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot_group3.txt'.format(halo), np.c_[group3_data[0], group3_data[1], group3_data[2], group3_data[3], group3_data[4], group3_data[5], group3_data[6], group3_data[7], #group3_data[8]], fmt='%1.16f', delimiter=' ', header='star_x_goodcenter kpc, star_y_goodcenter, star_r, star_smoothed_vrad km/s, star_smoothed_vrot, gas_vrad, gas_vrot,  star_dispersion')
#np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot_group4.txt'.format(halo), np.c_[group4_data[0], group4_data[1], group4_data[2], group4_data[3], group4_data[4], group4_data[5], group4_data[6], group4_data[7], #group4_data[8]], fmt='%1.16f', delimiter=' ', header='star_x_goodcenter kpc, star_y_goodcenter, star_r, star_smoothed_vrad km/s, star_smoothed_vrot, gas_vrad, gas_vrot,  star_dispersion')
#
#np.savetxt('/Volumes/FRIEND/analogs/data/{}_ad.txt'.format(halo), np.c_[group1_ad_error[0], group1_ad_error[1], group1_ad_error[2], group2_ad_error[0], group2_ad_error[1], group2_ad_error[2], group3_ad_error[0], group3_ad_error[1], group3_ad_error[2], group4_ad_error[0], group4_ad_error[1], group4_ad_error[2]], fmt='%1.16f', delimiter=' ', header='media, lower_error, upper_error; for group1, group2, group3, group4 in that order')

# ENIA: Are you gonna even be reading at this part of the page? Of course you will cos you are a tenacious, intelligent young woman that knows her code!! 
# ENIA: Don't be afraid to start over. It's a new chance to rebuild what you want <3 


# star1_ad = va(gas_avg_vrot, star1_avg_vrot)
# star2_ad = va(gas_avg_vrot, star2_avg_vrot)
# star3_ad = va(gas_avg_vrot, star3_avg_vrot)
# star4_ad = va(gas_avg_vrot, star4_avg_vrot)

# #remove nan values (that resulted from having empty radial bins)
# star1_ad = star1_ad[~np.isnan(star1_ad)]
# star2_ad = star2_ad[~np.isnan(star2_ad)]
# star3_ad = star3_ad[~np.isnan(star3_ad)]
# star4_ad = star4_ad[~np.isnan(star4_ad)]

#divide stars and gas into radial and average the rotation velocities within each bin -- want to change this to taking the median of AD in each bin
# R_min = 0 #kpc
# R_max = 21
# delta_r= 0.25 #kpc

# r_bins=np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)

# #will contain the average rotation velocity at each radial bin
# star1_avg_vrot = np.zeros(len(r_bins))
# star2_avg_vrot = np.zeros(len(r_bins))
# star3_avg_vrot = np.zeros(len(r_bins))
# star4_avg_vrot = np.zeros(len(r_bins))
# gas_avg_vrot = np.zeros(len(r_bins))

# #dividing each age groups into radial bin
# for i in range(len(r_bins)-1):
# 	star1_vrots=[b for b, a in zip(star1_vrot, star1_r) if a>=r_bins[i] and a<r_bins[i+1]]
# 	star1_avg_vrot[i]=np.median(star1_vrots)
# 	star2_vrots=[b for b, a in zip(star2_vrot, star2_r) if a>=r_bins[i] and a<r_bins[i+1]]
# 	star2_avg_vrot[i]=np.median(star2_vrots)
# 	star3_vrots=[b for b, a in zip(star3_vrot, star3_r) if a>=r_bins[i] and a<r_bins[i+1]]
# 	star3_avg_vrot[i]=np.median(star3_vrots)
# 	star4_vrots=[b for b, a in zip(star4_vrot, star4_r) if a>=r_bins[i] and a<r_bins[i+1]]
# 	star4_avg_vrot[i]=np.median(star4_vrots)
# 	gas_vrots=[b for b, a in zip(gas_vrot, gas_r) if a>=r_bins[i] and a<r_bins[i+1]]
# 	gas_avg_vrot[i]=np.median(gas_vrots)



# #plots-- 1 plot of rotation curves and 1 histogram of AD
# plt.scatter(r_bins[:-1], gas_avg_vrot[:-1], c = 'darkgrey', alpha = 0.4, label='gas')
# plt.scatter(r_bins[:-1], star1_avg_vrot[:-1], c = 'b', alpha = 0.4, label='30 Myr')
# plt.scatter(r_bins[:-1], star2_avg_vrot[:-1], c = 'm', alpha = 0.4, label='400 Myr')
# plt.scatter(r_bins[:-1], star3_avg_vrot[:-1], c = 'k', alpha = 0.4, label='2 Gyr')
# plt.scatter(r_bins[:-1], star4_avg_vrot[:-1], c = 'r', alpha = 0.4, label='4 Gyr')
# plt.ylim(0,300)
# plt.xlim(0,20)
# plt.ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=13)
# plt.xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
# plt.legend(loc=4)
# plt.savefig('/Volumes/FRIEND/analogs/plots/rc/{}_rotation_curve.pdf'.format(halo))
# plt.close()

# plt.hist(star1_ad[:-1], bins=range(-70, 70, 10), label='30 Myr', normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# plt.hist(star2_ad[:-1], bins=range(-70, 70, 10), label='400 Myr', normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m')
# plt.hist(star3_ad[:-1], bins=range(-70, 70, 10), label='2 Gyr', normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='k')
# plt.hist(star4_ad[:-1], bins=range(-70, 70, 10), label='4 Gyr', normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='r')
# plt.legend(loc=1, frameon=False)
# plt.xlim(-75,70)
# plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
# plt.savefig('/Volumes/FRIEND/analogs/plots/hist/{}_asymmetric_drift_hist.pdf'.format(halo))

# #save data to file
# np.savetxt('/Volumes/FRIEND/analogs/data/{}_vrot.txt'.format(halo), np.c_[r_bins[:-1], gas_avg_vrot[:-1], star1_avg_vrot[:-1], star2_avg_vrot[:-1], star3_avg_vrot[:-1], star4_avg_vrot[:-1]], fmt='%1.16f', delimiter=' ', header='r bin (kpc), avg v_rot gas (km/s), avg v_rot star 1 (km/s), avg v_rot star 2 (km/s), avg v_rot star 3 (km/s), avg v_rot star 4 (km/s)')

