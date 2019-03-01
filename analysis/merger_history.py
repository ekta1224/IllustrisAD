import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator

'''
reads in data from the halo properties files and from the files created in Illustris_curves.py
groups the analogs into bin based on their merger history
output: 2 plots (1 histogram of AD and one line graph of AD)
'''

halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt') #halo IDs
IDs, merger_times = np.loadtxt('../data/M31analogs_halo_props.txt', usecols=(0,9,), unpack=True) #halo properties file

#radial bins -- make sure below matches Illustris_curves.py ==============================================================================================
R_min = 0 #kpc
R_max = 21
delta_r = 0.1 #kpc
r_bins = np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)
r_bins = r_bins[:-1]
 
#count the number of mergers and non mergers =============================================================================================================
num_mergers = 0
num_no_mergers = 0

#will contain the number of nan values in each radial bin -- will be used to alter the averaging later
star1_num_nan_merger = np.zeros_like(r_bins)
star1_num_nan_no_merger = np.zeros_like(r_bins)
star2_num_nan_merger = np.zeros_like(r_bins)
star2_num_nan_no_merger = np.zeros_like(r_bins)
star3_num_nan_merger = np.zeros_like(r_bins)
star3_num_nan_no_merger = np.zeros_like(r_bins)
star4_num_nan_merger = np.zeros_like(r_bins)
star4_num_nan_no_merger = np.zeros_like(r_bins)
gas_num_nan_merger = np.zeros_like(r_bins)
gas_num_nan_no_merger = np.zeros_like(r_bins)

#will contain the summed vrots for each age bin and merger group
star1_vrot_merger = np.zeros_like(r_bins)
star1_vrot_no_merger = np.zeros_like(r_bins)
star2_vrot_merger = np.zeros_like(r_bins)
star2_vrot_no_merger = np.zeros_like(r_bins)
star3_vrot_merger = np.zeros_like(r_bins)
star3_vrot_no_merger = np.zeros_like(r_bins)
star4_vrot_merger = np.zeros_like(r_bins)
star4_vrot_no_merger = np.zeros_like(r_bins)
gas_vrot_merger = np.zeros_like(r_bins)
gas_vrot_no_merger = np.zeros_like(r_bins)
for halo in halos:
	#read in data
	print(halo)
	gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/FRIEND/analogs/data/{}_vrot.txt'.format(int(halo)), usecols=(1,2,3,4,5,), unpack=True) #km/s

	#get time of last merger
	N = np.where(halo == IDs)
	merger_time = merger_times[N]

	#divide into merger vs non merger groups
	if merger_time < 5: #Gyr
		#first populate num_nan so that the averaging is not affected by non entries
		for i in range(len(star1_vrot)):
			if np.isnan(star1_vrot[i]) == True:
				star1_num_nan_merger[i] = star1_num_nan_merger[i] + 1
			if np.isnan(star2_vrot[i]) == True:
				star2_num_nan_merger[i] = star2_num_nan_merger[i] + 1
			if np.isnan(star3_vrot[i]) == True:
				star3_num_nan_merger[i] = star3_num_nan_merger[i] + 1
			if np.isnan(star4_vrot[i]) == True:
				star4_num_nan_merger[i] = star4_num_nan_merger[i] + 1
			if np.isnan(gas_vrot[i]) == True:
				gas_num_nan_merger[i] = gas_num_nan_merger[i] + 1

		#convert nans to 0 so that array addition can happen
		star1_vrot = np.nan_to_num(star1_vrot)
		star2_vrot = np.nan_to_num(star2_vrot)
		star3_vrot = np.nan_to_num(star3_vrot)
		star4_vrot = np.nan_to_num(star4_vrot)
		gas_vrot = np.nan_to_num(gas_vrot)

		#add together all of the vrots per radial bin for all the halos
		star1_vrot_merger = star1_vrot_merger + star1_vrot
		star2_vrot_merger = star2_vrot_merger + star2_vrot
		star3_vrot_merger = star3_vrot_merger + star3_vrot
		star4_vrot_merger = star4_vrot_merger + star4_vrot
		gas_vrot_merger = gas_vrot_merger + gas_vrot
		num_mergers = num_mergers + 1

	else:
		#count how many nans there are so that the non data elements are not incorrectly averaged
		for i in range(len(star1_vrot)):
			if np.isnan(star1_vrot[i]) == True:
				star1_num_nan_no_merger[i] = star1_num_nan_no_merger[i] + 1
			if np.isnan(star2_vrot[i]) == True:
				star2_num_nan_no_merger[i] = star2_num_nan_no_merger[i] + 1
			if np.isnan(star3_vrot[i]) == True:
				star3_num_nan_no_merger[i] = star3_num_nan_no_merger[i] + 1
			if np.isnan(star4_vrot[i]) == True:
				star4_num_nan_no_merger[i] = star4_num_nan_no_merger[i] + 1
			if np.isnan(gas_vrot[i]) == True:
				gas_num_nan_no_merger[i] = gas_num_nan_no_merger[i] + 1

		#convert nans to 0 so that array addition can happen
		star1_vrot = np.nan_to_num(star1_vrot)
		star2_vrot = np.nan_to_num(star2_vrot)
		star3_vrot = np.nan_to_num(star3_vrot)
		star4_vrot = np.nan_to_num(star4_vrot)
		gas_vrot = np.nan_to_num(gas_vrot)

		#add together all of the vrots per radial bin for all the halos
		star1_vrot_no_merger = star1_vrot_no_merger + star1_vrot
		star2_vrot_no_merger = star2_vrot_no_merger + star2_vrot
		star3_vrot_no_merger = star3_vrot_no_merger + star3_vrot
		star4_vrot_no_merger = star4_vrot_no_merger + star4_vrot
		gas_vrot_no_merger = gas_vrot_no_merger + gas_vrot
		num_no_mergers = num_no_mergers + 1

#dividing the sums by the number of mergers - the number of nans in a given radial bin; now we have average rotation velocities ========
star1_vrot_merger = star1_vrot_merger / abs(star1_num_nan_merger - num_mergers)   
star1_vrot_no_merger = star1_vrot_no_merger / abs(star1_num_nan_no_merger - num_no_mergers)
star2_vrot_merger = star2_vrot_merger / abs(star2_num_nan_merger - num_mergers)     
star2_vrot_no_merger = star2_vrot_no_merger / abs(star2_num_nan_no_merger - num_no_mergers)
star3_vrot_merger = star3_vrot_merger / abs(star3_num_nan_merger - num_mergers)     
star3_vrot_no_merger = star3_vrot_no_merger / abs(star3_num_nan_no_merger - num_no_mergers)
star4_vrot_merger = star4_vrot_merger / abs(star4_num_nan_merger - num_mergers)     
star4_vrot_no_merger = star4_vrot_no_merger / abs(star4_num_nan_no_merger - num_no_mergers)
gas_vrot_merger = gas_vrot_merger / abs(gas_num_nan_merger - num_mergers)      
gas_vrot_no_merger = gas_vrot_no_merger / abs(gas_num_nan_no_merger - num_no_mergers) 

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

star1_AD_merger = va(gas_vrot_merger, star1_vrot_merger)
star2_AD_merger = va(gas_vrot_merger, star2_vrot_merger)
star3_AD_merger = va(gas_vrot_merger, star3_vrot_merger)
star4_AD_merger = va(gas_vrot_merger, star4_vrot_merger)
star1_AD_no_merger = va(gas_vrot_no_merger, star1_vrot_no_merger)
star2_AD_no_merger = va(gas_vrot_no_merger, star2_vrot_no_merger)
star3_AD_no_merger = va(gas_vrot_no_merger, star3_vrot_no_merger)
star4_AD_no_merger = va(gas_vrot_no_merger, star4_vrot_no_merger)

#get rid of nan values for the histograms
star1_AD_merger = star1_AD_merger[~np.isnan(star1_AD_merger)]
star2_AD_merger = star2_AD_merger[~np.isnan(star2_AD_merger)]
star3_AD_merger = star3_AD_merger[~np.isnan(star3_AD_merger)]
star4_AD_merger = star4_AD_merger[~np.isnan(star4_AD_merger)]
star1_AD_no_merger = star1_AD_no_merger[~np.isnan(star1_AD_no_merger)]
star2_AD_no_merger = star2_AD_no_merger[~np.isnan(star2_AD_no_merger)]
star3_AD_no_merger = star3_AD_no_merger[~np.isnan(star3_AD_no_merger)]
star4_AD_no_merger = star4_AD_no_merger[~np.isnan(star4_AD_no_merger)]

#plots ====================================================================================================================================
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

#AD histograms
f, axes= plt.subplots(1,2, sharey=True, sharex=False, figsize=(15,6))
for ax in axes:
	ax.set_xlim(-300, 300)
	ax.set_xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
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

axes[0].hist(star1_AD_merger, bins=range(-200, 300, 20), label='Group 1: Median AD = {}'.format(np.median(star1_AD_merger)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
axes[0].hist(star2_AD_merger, bins=range(-200, 300, 20), label='Group 2: Median AD = {}'.format(np.median(star2_AD_merger)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
axes[0].hist(star3_AD_merger, bins=range(-200, 300, 20), label='Group 3: Median AD = {}'.format(np.median(star3_AD_merger)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
axes[0].hist(star4_AD_merger, bins=range(-200, 300, 20), label='Group 4: Median AD = {}'.format(np.median(star4_AD_merger)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
axes[1].hist(star1_AD_no_merger, bins=range(-200, 300, 20), label='Group 1: Median AD = {}'.format(np.median(star1_AD_no_merger)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
axes[1].hist(star2_AD_no_merger, bins=range(-200, 300, 20), label='Group 2: Median AD = {}'.format(np.median(star2_AD_no_merger)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
axes[1].hist(star3_AD_no_merger, bins=range(-200, 300, 20), label='Group 3: Median AD = {}'.format(np.median(star3_AD_no_merger)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
axes[1].hist(star4_AD_no_merger, bins=range(-200, 300, 20), label='Group 4: Median AD = {}'.format(np.median(star4_AD_no_merger)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
axes[0].legend(loc=2, frameon=False)
axes[1].legend(loc=2, frameon=False)
axes[0].annotate('Recent 4:1 Merger', xy=(-270, 0.002), horizontalalignment='left', fontsize=12)
axes[1].annotate('No Recent 4:1 Merger', xy=(-270, 0.002), horizontalalignment='left', fontsize=12)
plt.subplots_adjust(wspace=.09, hspace=0)
plt.savefig('/Volumes/FRIEND/analogs/plots/merger/AD.png', bbox_inches='tight')
plt.close()

#RCs
f, axes= plt.subplots(1,2, sharey=True, sharex=False, figsize=(15,6))
for ax in axes:
	ax.set_xlim(0, 20)
	ax.set_xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
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
axes[0].set_ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=13)
ax.set_ylim(0, 500)

axes[0].scatter(r_bins, gas_vrot_merger, c = 'darkgrey', s=12, label='gas')
axes[0].scatter(r_bins, star1_vrot_merger, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1')
axes[0].scatter(r_bins, star2_vrot_merger, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2')
axes[0].scatter(r_bins, star3_vrot_merger, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3')
axes[0].scatter(r_bins, star4_vrot_merger, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4')
axes[1].scatter(r_bins, gas_vrot_no_merger, c = 'darkgrey', s=12, label='gas')
axes[1].scatter(r_bins, star1_vrot_no_merger, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1')
axes[1].scatter(r_bins, star2_vrot_no_merger, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2')
axes[1].scatter(r_bins, star3_vrot_no_merger, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3')
axes[1].scatter(r_bins, star4_vrot_no_merger, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4')
axes[0].legend(loc=2, frameon=False)
axes[0].annotate('Recent 4:1 Merger', xy=(12.5, 460), horizontalalignment='left', fontsize=12)
axes[1].annotate('No Recent 4:1 Merger', xy=(12.5, 460), horizontalalignment='left', fontsize=12)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Volumes/FRIEND/analogs/plots/merger/rc.png', bbox_inches='tight')
plt.close()

#line plot
#need the nan AD values so same dimension as the rbins
star1_AD_merger = va(gas_vrot_merger, star1_vrot_merger)
star2_AD_merger = va(gas_vrot_merger, star2_vrot_merger)
star3_AD_merger = va(gas_vrot_merger, star3_vrot_merger)
star4_AD_merger = va(gas_vrot_merger, star4_vrot_merger)
star1_AD_no_merger = va(gas_vrot_no_merger, star1_vrot_no_merger)
star2_AD_no_merger = va(gas_vrot_no_merger, star2_vrot_no_merger)
star3_AD_no_merger = va(gas_vrot_no_merger, star3_vrot_no_merger)
star4_AD_no_merger = va(gas_vrot_no_merger, star4_vrot_no_merger)

f, axes= plt.subplots(1,2, sharey=True, sharex=False, figsize=(15,6))
for ax in axes:
	ax.set_xlim(0, 20)
	ax.set_xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
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
axes[0].set_ylabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
ax.set_ylim(-300, 300)

axes[0].scatter(r_bins, star1_AD_merger, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1')
axes[0].scatter(r_bins, star2_AD_merger, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2')
axes[0].scatter(r_bins, star3_AD_merger, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3')
axes[0].scatter(r_bins, star4_AD_merger, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4')
axes[1].scatter(r_bins, star1_AD_no_merger, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1')
axes[1].scatter(r_bins, star2_AD_no_merger, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2')
axes[1].scatter(r_bins, star3_AD_no_merger, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3')
axes[1].scatter(r_bins, star4_AD_no_merger, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4')
axes[0].legend(loc=2, frameon=False)
axes[0].annotate('Recent 4:1 Merger', xy=(12.5, -260), horizontalalignment='left', fontsize=12)
axes[1].annotate('No Recent 4:1 Merger', xy=(12.5, -260), horizontalalignment='left', fontsize=12)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Volumes/FRIEND/analogs/plots/merger/ad_lineplot.png', bbox_inches='tight')
plt.close()







