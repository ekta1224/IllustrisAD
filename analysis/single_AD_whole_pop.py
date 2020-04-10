import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc
'''
reads in data from the halo properties files and from the files created in Illustris_curves.py
groups the analogs into bin based on their merger history
output: 2 plots (1 histogram of AD and one line graph of AD)
'''

halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt') #halo IDs M31analog_IDs_IllustrisAD.txt
IDs, merger_times = np.loadtxt('../data/M31analogs_halo_props.txt', usecols=(0,10,), unpack=True) #halo properties file
# noM33_ids = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/M31_analogs_IDs_noM33_TNG100.txt')
# M33_ids = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/M31_analogs_IDs_M33_TNG100.txt')
# halos = [int(a) for a in M33_ids]



#radial bins -- make sure below matches Illustris_curves.py ==============================================================================================
R_min = 2 #kpc #CHANGE THIS TO 2 WHEN DOING SMOOTHED
R_max = 21
delta_r = 0.1 #kpc
r_bins = np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)
r_bins = r_bins[:-1]
 
#count the number of mergers and non mergers =============================================================================================================
num_mergers = 0
num_no_mergers = 0

#contains non nan data to be used for the average and the standard deviation
star1_data_merger = [[] for i in range(len(r_bins))]
star2_data_merger = [[] for i in range(len(r_bins))]
star3_data_merger = [[] for i in range(len(r_bins))]
star4_data_merger = [[] for i in range(len(r_bins))]
gas_data_merger = [[] for i in range(len(r_bins))]
star1_data_no_merger = [[] for i in range(len(r_bins))]
star2_data_no_merger = [[] for i in range(len(r_bins))]
star3_data_no_merger = [[] for i in range(len(r_bins))]
star4_data_no_merger = [[] for i in range(len(r_bins))]
gas_data_no_merger = [[] for i in range(len(r_bins))]
for halo in halos:
	#read in data
	print(halo)
	gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/Titan/analogs/data/{}_vrot_gas_smooted.txt'.format(int(halo)), usecols=(1,2,3,4,5,), unpack=True) #km/s

	#get time of last merger
	N = np.where(halo == IDs)
	merger_time = merger_times[N]

	#divide into merger vs non merger groups
	if merger_time <= 14:#np.median(merger_times): 
		for i in range(len(star1_vrot)):
			if np.isnan(star1_vrot[i]) == False:
				star1_data_merger[i].append(star1_vrot[i])
			if np.isnan(star2_vrot[i]) == False:
				star2_data_merger[i].append(star2_vrot[i])
			if np.isnan(star3_vrot[i]) == False:
				star3_data_merger[i].append(star3_vrot[i])
			if np.isnan(star4_vrot[i]) == False:
				star4_data_merger[i].append(star4_vrot[i])
			if np.isnan(gas_vrot[i]) == False:
				gas_data_merger[i].append(gas_vrot[i])
		num_mergers = num_mergers + 1

	else:
		#count how many nans there are so that the non data elements are not incorrectly averaged
		for i in range(len(star1_vrot)):
			if np.isnan(star1_vrot[i]) == False:
				star1_data_no_merger[i].append(star1_vrot[i])
			if np.isnan(star2_vrot[i]) == False:
				star2_data_no_merger[i].append(star2_vrot[i])
			if np.isnan(star3_vrot[i]) == False:
				star3_data_no_merger[i].append(star3_vrot[i])
			if np.isnan(star4_vrot[i]) == False:
				star4_data_no_merger[i].append(star4_vrot[i])
			if np.isnan(gas_vrot[i]) == False:
				gas_data_no_merger[i].append(gas_vrot[i])
		num_no_mergers = num_no_mergers + 1

#averages
def calc_errors(data):
	if len(data) == 0:
		return np.nan, np.nan, np.nan
	else:
		result = np.percentile(data, [16, 50, 84])
		median = result[1]
		lower_error = result[1] - result[0]
		upper_error = result[2] - result[1]
		return median, lower_error, upper_error


star1_vrot_merger = np.zeros_like(r_bins)
star2_vrot_merger = np.zeros_like(r_bins)
star3_vrot_merger = np.zeros_like(r_bins)
star4_vrot_merger = np.zeros_like(r_bins)
gas_vrot_merger = np.zeros_like(r_bins)
star1_vrot_no_merger = np.zeros_like(r_bins)
star2_vrot_no_merger = np.zeros_like(r_bins)
star3_vrot_no_merger = np.zeros_like(r_bins)
star4_vrot_no_merger = np.zeros_like(r_bins)
gas_vrot_no_merger = np.zeros_like(r_bins)

star1_upper_error_merger = np.zeros_like(r_bins)
star2_upper_error_merger = np.zeros_like(r_bins)
star3_upper_error_merger = np.zeros_like(r_bins)
star4_upper_error_merger = np.zeros_like(r_bins)
gas_upper_error_merger = np.zeros_like(r_bins)
star1_upper_error_no_merger = np.zeros_like(r_bins)
star2_upper_error_no_merger = np.zeros_like(r_bins)
star3_upper_error_no_merger = np.zeros_like(r_bins)
star4_upper_error_no_merger = np.zeros_like(r_bins)
gas_upper_error_no_merger = np.zeros_like(r_bins)

star1_lower_error_merger = np.zeros_like(r_bins)
star2_lower_error_merger = np.zeros_like(r_bins)
star3_lower_error_merger = np.zeros_like(r_bins)
star4_lower_error_merger = np.zeros_like(r_bins)
gas_lower_error_merger = np.zeros_like(r_bins)
star1_lower_error_no_merger = np.zeros_like(r_bins)
star2_lower_error_no_merger = np.zeros_like(r_bins)
star3_lower_error_no_merger = np.zeros_like(r_bins)
star4_lower_error_no_merger = np.zeros_like(r_bins)
gas_lower_error_no_merger = np.zeros_like(r_bins)

#median and standard deviation
for i in range(len(r_bins)):
	#errors on the ROTATION velocity -- is this same for AD?
	star1_vrot_merger[i], star1_lower_error_merger[i], star1_upper_error_merger[i] = calc_errors(star1_data_merger[i])
	star2_vrot_merger[i], star2_lower_error_merger[i], star2_upper_error_merger[i] = calc_errors(star2_data_merger[i])
	star3_vrot_merger[i], star3_lower_error_merger[i], star3_upper_error_merger[i] = calc_errors(star3_data_merger[i])
	star4_vrot_merger[i], star4_lower_error_merger[i], star4_upper_error_merger[i] = calc_errors(star4_data_merger[i])
	gas_vrot_merger[i], gas_lower_error_merger[i], gas_upper_error_merger[i] = calc_errors(gas_data_merger[i])
	star1_vrot_no_merger[i], star1_lower_error_no_merger[i], star1_upper_error_no_merger[i] = calc_errors(star1_data_no_merger[i])
	star2_vrot_no_merger[i], star2_lower_error_no_merger[i], star2_upper_error_no_merger[i] = calc_errors(star2_data_no_merger[i])
	star3_vrot_no_merger[i], star3_lower_error_no_merger[i], star3_upper_error_no_merger[i] = calc_errors(star3_data_no_merger[i])
	star4_vrot_no_merger[i], star4_lower_error_no_merger[i], star4_upper_error_no_merger[i] = calc_errors(star4_data_no_merger[i])
	gas_vrot_no_merger[i], gas_lower_error_no_merger[i], gas_upper_error_no_merger[i] = calc_errors(gas_data_no_merger[i])

#should I divide this by the number in each group?
star1_error_merger = np.row_stack((star1_lower_error_merger, star1_upper_error_merger))
star2_error_merger = np.row_stack((star2_lower_error_merger, star2_upper_error_merger))
star3_error_merger = np.row_stack((star3_lower_error_merger, star3_upper_error_merger))
star4_error_merger = np.row_stack((star4_lower_error_merger, star4_upper_error_merger))
gas_error_merger = np.row_stack((gas_lower_error_merger, gas_upper_error_merger)) 

star1_error_no_merger = np.row_stack((star1_lower_error_no_merger, star1_upper_error_no_merger))
star2_error_no_merger = np.row_stack((star2_lower_error_no_merger, star2_upper_error_no_merger))
star3_error_no_merger = np.row_stack((star3_lower_error_no_merger, star3_upper_error_no_merger))
star4_error_no_merger = np.row_stack((star4_lower_error_no_merger, star4_upper_error_no_merger))
gas_error_no_merger = np.row_stack((gas_lower_error_no_merger, gas_upper_error_no_merger))

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
	        ax.spines[axis].set_linewidth(1)
	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=1)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()

#AD histograms
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

axes[1].hist(star1_AD_merger, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_merger),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
axes[1].hist(star2_AD_merger, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_merger),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
axes[1].hist(star3_AD_merger, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_merger),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
axes[1].hist(star4_AD_merger, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_merger),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
axes[1].legend(loc=1, frameon=False, fontsize=12)
axes[1].set_xlim(-115, 150)
axes[1].set_xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=16)
axes[1].set_ylabel(r'$ \rm PDF$', fontsize=16)
#single_plot()
axes[0].errorbar(r_bins, gas_vrot_merger,   yerr = gas_error_merger, errorevery=20, c = 'darkgrey', label='gas')
axes[0].errorbar(r_bins, star1_vrot_merger, yerr = star1_error_merger, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
axes[0].errorbar(r_bins, star2_vrot_merger, yerr = star2_error_merger, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
axes[0].errorbar(r_bins, star3_vrot_merger, yerr = star3_error_merger, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
axes[0].errorbar(r_bins, star4_vrot_merger, yerr = star4_error_merger, errorevery=17, c = 'r', alpha = 0.6, linestyle=':', label=r'$\rm \geq 10\ Gyr$')
axes[0].set_ylim(0,220)
axes[0].set_xlim(0,22)
axes[0].set_ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=16)
axes[0].set_xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=16)
axes[0].legend(loc=2, frameon=False, fontsize=11)

#plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('/Users/amandaquirk/Desktop/wholepop_AD_RC.png', bbox_inches='tight')
plt.close()
