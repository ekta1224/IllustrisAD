import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc 
from matplotlib.ticker import MaxNLocator
'''
reads in data from the halo properties files and from the files created in Illustris_curves.py
groups the analogs into bin based on their merger history
output: 2 plots (1 histogram of AD and one line graph of AD)
'''

halos = np.loadtxt('/Volumes/Titan/analogs/IllustrisAD/TNG_v2/M31analogs_halo_props_TNG100_revised.txt', usecols=(0,), unpack=True) #halo IDs M31analog_IDs_IllustrisAD.txt
IDs, merger_times = np.loadtxt('/Volumes/Titan/analogs/IllustrisAD/TNG_v2/M31analogs_merger_props_TNG100_revised.txt', usecols=(0,2,), unpack=True) #halo properties file
# noM33_ids = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/M31_analogs_IDs_noM33_TNG100.txt')
# M33_ids = np.loadtxt('/Volumes/FRIEND/analogs/TNGdata/M31_analogs_IDs_M33_TNG100.txt')
# halos = [int(a) for a in M33_ids]
med_time = 7.99
merger_time_division = [4, med_time, 12]

#radial bins -- make sure below matches Illustris_curves.py ==============================================================================================
R_min = 2 #kpc #CHANGE THIS TO 2 WHEN DOING SMOOTHED
R_max = 21
delta_r = 0.1 #kpc
r_bins = np.linspace(R_min, R_max, (R_max - R_min) / delta_r + 1)
r_bins = r_bins[:-1]
 
#count the number of mergers and non mergers =============================================================================================================
num_mergers1 = 0
num_no_mergers1 = 0

num_mergers2 = 0
num_no_mergers2 = 0

num_mergers3 = 0
num_no_mergers3 = 0
#contains non nan data to be used for the average and the standard deviation
star1_data_merger1 = [[] for i in range(len(r_bins))]
star2_data_merger1 = [[] for i in range(len(r_bins))]
star3_data_merger1 = [[] for i in range(len(r_bins))]
star4_data_merger1 = [[] for i in range(len(r_bins))]
gas_data_merger1 = [[] for i in range(len(r_bins))]
star1_data_no_merger1 = [[] for i in range(len(r_bins))]
star2_data_no_merger1 = [[] for i in range(len(r_bins))]
star3_data_no_merger1 = [[] for i in range(len(r_bins))]
star4_data_no_merger1 = [[] for i in range(len(r_bins))]
gas_data_no_merger1 = [[] for i in range(len(r_bins))]

star1_data_merger2 = [[] for i in range(len(r_bins))]
star2_data_merger2 = [[] for i in range(len(r_bins))]
star3_data_merger2 = [[] for i in range(len(r_bins))]
star4_data_merger2 = [[] for i in range(len(r_bins))]
gas_data_merger2 = [[] for i in range(len(r_bins))]
star1_data_no_merger2 = [[] for i in range(len(r_bins))]
star2_data_no_merger2 = [[] for i in range(len(r_bins))]
star3_data_no_merger2 = [[] for i in range(len(r_bins))]
star4_data_no_merger2 = [[] for i in range(len(r_bins))]
gas_data_no_merger2 = [[] for i in range(len(r_bins))]

star1_data_merger3 = [[] for i in range(len(r_bins))]
star2_data_merger3 = [[] for i in range(len(r_bins))]
star3_data_merger3 = [[] for i in range(len(r_bins))]
star4_data_merger3 = [[] for i in range(len(r_bins))]
gas_data_merger3 = [[] for i in range(len(r_bins))]
star1_data_no_merger3 = [[] for i in range(len(r_bins))]
star2_data_no_merger3 = [[] for i in range(len(r_bins))]
star3_data_no_merger3 = [[] for i in range(len(r_bins))]
star4_data_no_merger3 = [[] for i in range(len(r_bins))]
gas_data_no_merger3 = [[] for i in range(len(r_bins))]
for halo in halos:
	#read in data
	print(halo)
	gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/Titan/analogs/TNGdata/smoothed_vrot/{}_vrot_smoothed.txt'.format(int(halo)), usecols=(1,2,3,4,5,), unpack=True) #km/s

	#get time of last merger
	N = np.where(halo == IDs)
	merger_time = merger_times[N]

	#divide into merger vs non merger groups
	for time in merger_time_division:
		if time == 4:
			if merger_time <= time:#np.median(merger_times): 
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_merger1[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_merger1[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_merger1[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_merger1[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_merger1[i].append(gas_vrot[i])
				num_mergers1 = num_mergers1 + 1
		
			else:
				#count how many nans there are so that the non data elements are not incorrectly averaged
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_no_merger1[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_no_merger1[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_no_merger1[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_no_merger1[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_no_merger1[i].append(gas_vrot[i])
				num_no_mergers1 = num_no_mergers1 + 1
		elif time == med_time:
			if merger_time <= time:#np.median(merger_times): 
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_merger2[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_merger2[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_merger2[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_merger2[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_merger2[i].append(gas_vrot[i])
				num_mergers2 = num_mergers2 + 1
		
			else:
				#count how many nans there are so that the non data elements are not incorrectly averaged
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_no_merger2[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_no_merger2[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_no_merger2[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_no_merger2[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_no_merger2[i].append(gas_vrot[i])
				num_no_mergers2 = num_no_mergers2 + 1
		else:
			if merger_time <= time:#np.median(merger_times): 
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_merger3[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_merger3[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_merger3[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_merger3[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_merger3[i].append(gas_vrot[i])
				num_mergers3 = num_mergers3 + 1
		
			else:
				#count how many nans there are so that the non data elements are not incorrectly averaged
				for i in range(len(star1_vrot)):
					if np.isnan(star1_vrot[i]) == False:
						star1_data_no_merger3[i].append(star1_vrot[i])
					if np.isnan(star2_vrot[i]) == False:
						star2_data_no_merger3[i].append(star2_vrot[i])
					if np.isnan(star3_vrot[i]) == False:
						star3_data_no_merger3[i].append(star3_vrot[i])
					if np.isnan(star4_vrot[i]) == False:
						star4_data_no_merger3[i].append(star4_vrot[i])
					if np.isnan(gas_vrot[i]) == False:
						gas_data_no_merger3[i].append(gas_vrot[i])
				num_no_mergers3 = num_no_mergers3 + 1

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


star1_vrot_merger1 = np.zeros_like(r_bins)
star2_vrot_merger1 = np.zeros_like(r_bins)
star3_vrot_merger1 = np.zeros_like(r_bins)
star4_vrot_merger1 = np.zeros_like(r_bins)
gas_vrot_merger1 = np.zeros_like(r_bins)
star1_vrot_no_merger1 = np.zeros_like(r_bins)
star2_vrot_no_merger1 = np.zeros_like(r_bins)
star3_vrot_no_merger1 = np.zeros_like(r_bins)
star4_vrot_no_merger1 = np.zeros_like(r_bins)
gas_vrot_no_merger1 = np.zeros_like(r_bins)

star1_upper_error_merger1 = np.zeros_like(r_bins)
star2_upper_error_merger1 = np.zeros_like(r_bins)
star3_upper_error_merger1 = np.zeros_like(r_bins)
star4_upper_error_merger1 = np.zeros_like(r_bins)
gas_upper_error_merger1 = np.zeros_like(r_bins)
star1_upper_error_no_merger1 = np.zeros_like(r_bins)
star2_upper_error_no_merger1 = np.zeros_like(r_bins)
star3_upper_error_no_merger1 = np.zeros_like(r_bins)
star4_upper_error_no_merger1 = np.zeros_like(r_bins)
gas_upper_error_no_merger1 = np.zeros_like(r_bins)

star1_lower_error_merger1 = np.zeros_like(r_bins)
star2_lower_error_merger1 = np.zeros_like(r_bins)
star3_lower_error_merger1 = np.zeros_like(r_bins)
star4_lower_error_merger1 = np.zeros_like(r_bins)
gas_lower_error_merger1 = np.zeros_like(r_bins)
star1_lower_error_no_merger1 = np.zeros_like(r_bins)
star2_lower_error_no_merger1 = np.zeros_like(r_bins)
star3_lower_error_no_merger1 = np.zeros_like(r_bins)
star4_lower_error_no_merger1 = np.zeros_like(r_bins)
gas_lower_error_no_merger1 = np.zeros_like(r_bins)

star1_vrot_merger2 = np.zeros_like(r_bins)
star2_vrot_merger2 = np.zeros_like(r_bins)
star3_vrot_merger2 = np.zeros_like(r_bins)
star4_vrot_merger2 = np.zeros_like(r_bins)
gas_vrot_merger2 = np.zeros_like(r_bins)
star1_vrot_no_merger2 = np.zeros_like(r_bins)
star2_vrot_no_merger2 = np.zeros_like(r_bins)
star3_vrot_no_merger2 = np.zeros_like(r_bins)
star4_vrot_no_merger2 = np.zeros_like(r_bins)
gas_vrot_no_merger2 = np.zeros_like(r_bins)

star1_upper_error_merger2 = np.zeros_like(r_bins)
star2_upper_error_merger2 = np.zeros_like(r_bins)
star3_upper_error_merger2 = np.zeros_like(r_bins)
star4_upper_error_merger2 = np.zeros_like(r_bins)
gas_upper_error_merger2 = np.zeros_like(r_bins)
star1_upper_error_no_merger2 = np.zeros_like(r_bins)
star2_upper_error_no_merger2 = np.zeros_like(r_bins)
star3_upper_error_no_merger2 = np.zeros_like(r_bins)
star4_upper_error_no_merger2 = np.zeros_like(r_bins)
gas_upper_error_no_merger2 = np.zeros_like(r_bins)

star1_lower_error_merger2 = np.zeros_like(r_bins)
star2_lower_error_merger2 = np.zeros_like(r_bins)
star3_lower_error_merger2 = np.zeros_like(r_bins)
star4_lower_error_merger2 = np.zeros_like(r_bins)
gas_lower_error_merger2 = np.zeros_like(r_bins)
star1_lower_error_no_merger2 = np.zeros_like(r_bins)
star2_lower_error_no_merger2 = np.zeros_like(r_bins)
star3_lower_error_no_merger2 = np.zeros_like(r_bins)
star4_lower_error_no_merger2 = np.zeros_like(r_bins)
gas_lower_error_no_merger2 = np.zeros_like(r_bins)

star1_vrot_merger3 = np.zeros_like(r_bins)
star2_vrot_merger3 = np.zeros_like(r_bins)
star3_vrot_merger3 = np.zeros_like(r_bins)
star4_vrot_merger3 = np.zeros_like(r_bins)
gas_vrot_merger3 = np.zeros_like(r_bins)
star1_vrot_no_merger3 = np.zeros_like(r_bins)
star2_vrot_no_merger3 = np.zeros_like(r_bins)
star3_vrot_no_merger3 = np.zeros_like(r_bins)
star4_vrot_no_merger3 = np.zeros_like(r_bins)
gas_vrot_no_merger3 = np.zeros_like(r_bins)

star1_upper_error_merger3 = np.zeros_like(r_bins)
star2_upper_error_merger3 = np.zeros_like(r_bins)
star3_upper_error_merger3 = np.zeros_like(r_bins)
star4_upper_error_merger3 = np.zeros_like(r_bins)
gas_upper_error_merger3 = np.zeros_like(r_bins)
star1_upper_error_no_merger3 = np.zeros_like(r_bins)
star2_upper_error_no_merger3 = np.zeros_like(r_bins)
star3_upper_error_no_merger3 = np.zeros_like(r_bins)
star4_upper_error_no_merger3 = np.zeros_like(r_bins)
gas_upper_error_no_merger3 = np.zeros_like(r_bins)

star1_lower_error_merger3 = np.zeros_like(r_bins)
star2_lower_error_merger3 = np.zeros_like(r_bins)
star3_lower_error_merger3 = np.zeros_like(r_bins)
star4_lower_error_merger3 = np.zeros_like(r_bins)
gas_lower_error_merger3 = np.zeros_like(r_bins)
star1_lower_error_no_merger3 = np.zeros_like(r_bins)
star2_lower_error_no_merger3 = np.zeros_like(r_bins)
star3_lower_error_no_merger3 = np.zeros_like(r_bins)
star4_lower_error_no_merger3 = np.zeros_like(r_bins)
gas_lower_error_no_merger3 = np.zeros_like(r_bins)

#median and standard deviation
for i in range(len(r_bins)):
	#errors on the ROTATION velocity -- is this same for AD?
	star1_vrot_merger1[i], star1_lower_error_merger1[i], star1_upper_error_merger1[i] = calc_errors(star1_data_merger1[i])
	star2_vrot_merger1[i], star2_lower_error_merger1[i], star2_upper_error_merger1[i] = calc_errors(star2_data_merger1[i])
	star3_vrot_merger1[i], star3_lower_error_merger1[i], star3_upper_error_merger1[i] = calc_errors(star3_data_merger1[i])
	star4_vrot_merger1[i], star4_lower_error_merger1[i], star4_upper_error_merger1[i] = calc_errors(star4_data_merger1[i])
	gas_vrot_merger1[i], gas_lower_error_merger1[i], gas_upper_error_merger1[i] = calc_errors(gas_data_merger1[i])
	star1_vrot_no_merger1[i], star1_lower_error_no_merger1[i], star1_upper_error_no_merger1[i] = calc_errors(star1_data_no_merger1[i])
	star2_vrot_no_merger1[i], star2_lower_error_no_merger1[i], star2_upper_error_no_merger1[i] = calc_errors(star2_data_no_merger1[i])
	star3_vrot_no_merger1[i], star3_lower_error_no_merger1[i], star3_upper_error_no_merger1[i] = calc_errors(star3_data_no_merger1[i])
	star4_vrot_no_merger1[i], star4_lower_error_no_merger1[i], star4_upper_error_no_merger1[i] = calc_errors(star4_data_no_merger1[i])
	gas_vrot_no_merger1[i], gas_lower_error_no_merger1[i], gas_upper_error_no_merger1[i] = calc_errors(gas_data_no_merger1[i])

	star1_vrot_merger2[i], star1_lower_error_merger2[i], star1_upper_error_merger2[i] = calc_errors(star1_data_merger2[i])
	star2_vrot_merger2[i], star2_lower_error_merger2[i], star2_upper_error_merger2[i] = calc_errors(star2_data_merger2[i])
	star3_vrot_merger2[i], star3_lower_error_merger2[i], star3_upper_error_merger2[i] = calc_errors(star3_data_merger2[i])
	star4_vrot_merger2[i], star4_lower_error_merger2[i], star4_upper_error_merger2[i] = calc_errors(star4_data_merger2[i])
	gas_vrot_merger2[i], gas_lower_error_merger2[i], gas_upper_error_merger2[i] = calc_errors(gas_data_merger2[i])
	star1_vrot_no_merger2[i], star1_lower_error_no_merger2[i], star1_upper_error_no_merger2[i] = calc_errors(star1_data_no_merger2[i])
	star2_vrot_no_merger2[i], star2_lower_error_no_merger2[i], star2_upper_error_no_merger2[i] = calc_errors(star2_data_no_merger2[i])
	star3_vrot_no_merger2[i], star3_lower_error_no_merger2[i], star3_upper_error_no_merger2[i] = calc_errors(star3_data_no_merger2[i])
	star4_vrot_no_merger2[i], star4_lower_error_no_merger2[i], star4_upper_error_no_merger2[i] = calc_errors(star4_data_no_merger2[i])
	gas_vrot_no_merger2[i], gas_lower_error_no_merger2[i], gas_upper_error_no_merger2[i] = calc_errors(gas_data_no_merger2[i])

	star1_vrot_merger3[i], star1_lower_error_merger3[i], star1_upper_error_merger3[i] = calc_errors(star1_data_merger3[i])
	star2_vrot_merger3[i], star2_lower_error_merger3[i], star2_upper_error_merger3[i] = calc_errors(star2_data_merger3[i])
	star3_vrot_merger3[i], star3_lower_error_merger3[i], star3_upper_error_merger3[i] = calc_errors(star3_data_merger3[i])
	star4_vrot_merger3[i], star4_lower_error_merger3[i], star4_upper_error_merger3[i] = calc_errors(star4_data_merger3[i])
	gas_vrot_merger3[i], gas_lower_error_merger3[i], gas_upper_error_merger3[i] = calc_errors(gas_data_merger3[i])
	star1_vrot_no_merger3[i], star1_lower_error_no_merger3[i], star1_upper_error_no_merger3[i] = calc_errors(star1_data_no_merger3[i])
	star2_vrot_no_merger3[i], star2_lower_error_no_merger3[i], star2_upper_error_no_merger3[i] = calc_errors(star2_data_no_merger3[i])
	star3_vrot_no_merger3[i], star3_lower_error_no_merger3[i], star3_upper_error_no_merger3[i] = calc_errors(star3_data_no_merger3[i])
	star4_vrot_no_merger3[i], star4_lower_error_no_merger3[i], star4_upper_error_no_merger3[i] = calc_errors(star4_data_no_merger3[i])
	gas_vrot_no_merger3[i], gas_lower_error_no_merger3[i], gas_upper_error_no_merger3[i] = calc_errors(gas_data_no_merger3[i])

star1_error_merger1 = np.row_stack((star1_lower_error_merger1, star1_upper_error_merger1))
star2_error_merger1 = np.row_stack((star2_lower_error_merger1, star2_upper_error_merger1))
star3_error_merger1 = np.row_stack((star3_lower_error_merger1, star3_upper_error_merger1))
star4_error_merger1 = np.row_stack((star4_lower_error_merger1, star4_upper_error_merger1))
gas_error_merger1 = np.row_stack((gas_lower_error_merger1, gas_upper_error_merger1)) 

star1_error_no_merger1 = np.row_stack((star1_lower_error_no_merger1, star1_upper_error_no_merger1))
star2_error_no_merger1 = np.row_stack((star2_lower_error_no_merger1, star2_upper_error_no_merger1))
star3_error_no_merger1 = np.row_stack((star3_lower_error_no_merger1, star3_upper_error_no_merger1))
star4_error_no_merger1 = np.row_stack((star4_lower_error_no_merger1, star4_upper_error_no_merger1))
gas_error_no_merger1 = np.row_stack((gas_lower_error_no_merger1, gas_upper_error_no_merger1))

star1_error_merger2 = np.row_stack((star1_lower_error_merger2, star1_upper_error_merger2))
star2_error_merger2 = np.row_stack((star2_lower_error_merger2, star2_upper_error_merger2))
star3_error_merger2 = np.row_stack((star3_lower_error_merger2, star3_upper_error_merger2))
star4_error_merger2 = np.row_stack((star4_lower_error_merger2, star4_upper_error_merger2))
gas_error_merger2 = np.row_stack((gas_lower_error_merger2, gas_upper_error_merger2)) 

star1_error_no_merger2 = np.row_stack((star1_lower_error_no_merger2, star1_upper_error_no_merger2))
star2_error_no_merger2 = np.row_stack((star2_lower_error_no_merger2, star2_upper_error_no_merger2))
star3_error_no_merger2 = np.row_stack((star3_lower_error_no_merger2, star3_upper_error_no_merger2))
star4_error_no_merger2 = np.row_stack((star4_lower_error_no_merger2, star4_upper_error_no_merger2))
gas_error_no_merger2 = np.row_stack((gas_lower_error_no_merger2, gas_upper_error_no_merger2))

star1_error_merger3 = np.row_stack((star1_lower_error_merger3, star1_upper_error_merger3))
star2_error_merger3 = np.row_stack((star2_lower_error_merger3, star2_upper_error_merger3))
star3_error_merger3 = np.row_stack((star3_lower_error_merger3, star3_upper_error_merger3))
star4_error_merger3 = np.row_stack((star4_lower_error_merger3, star4_upper_error_merger3))
gas_error_merger3 = np.row_stack((gas_lower_error_merger3, gas_upper_error_merger3)) 

star1_error_no_merger3 = np.row_stack((star1_lower_error_no_merger3, star1_upper_error_no_merger3))
star2_error_no_merger3 = np.row_stack((star2_lower_error_no_merger3, star2_upper_error_no_merger3))
star3_error_no_merger3 = np.row_stack((star3_lower_error_no_merger3, star3_upper_error_no_merger3))
star4_error_no_merger3 = np.row_stack((star4_lower_error_no_merger3, star4_upper_error_no_merger3))
gas_error_no_merger3 = np.row_stack((gas_lower_error_no_merger3, gas_upper_error_no_merger3))

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

star1_AD_merger1 = va(gas_vrot_merger1, star1_vrot_merger1)
star2_AD_merger1 = va(gas_vrot_merger1, star2_vrot_merger1)
star3_AD_merger1 = va(gas_vrot_merger1, star3_vrot_merger1)
star4_AD_merger1 = va(gas_vrot_merger1, star4_vrot_merger1)
star1_AD_no_merger1 = va(gas_vrot_no_merger1, star1_vrot_no_merger1)
star2_AD_no_merger1 = va(gas_vrot_no_merger1, star2_vrot_no_merger1)
star3_AD_no_merger1 = va(gas_vrot_no_merger1, star3_vrot_no_merger1)
star4_AD_no_merger1 = va(gas_vrot_no_merger1, star4_vrot_no_merger1)

star1_AD_merger2 = va(gas_vrot_merger2, star1_vrot_merger2)
star2_AD_merger2 = va(gas_vrot_merger2, star2_vrot_merger2)
star3_AD_merger2 = va(gas_vrot_merger2, star3_vrot_merger2)
star4_AD_merger2 = va(gas_vrot_merger2, star4_vrot_merger2)
star1_AD_no_merger2 = va(gas_vrot_no_merger2, star1_vrot_no_merger2)
star2_AD_no_merger2 = va(gas_vrot_no_merger2, star2_vrot_no_merger2)
star3_AD_no_merger2 = va(gas_vrot_no_merger2, star3_vrot_no_merger2)
star4_AD_no_merger2 = va(gas_vrot_no_merger2, star4_vrot_no_merger2)

star1_AD_merger3 = va(gas_vrot_merger3, star1_vrot_merger3)
star2_AD_merger3 = va(gas_vrot_merger3, star2_vrot_merger3)
star3_AD_merger3 = va(gas_vrot_merger3, star3_vrot_merger3)
star4_AD_merger3 = va(gas_vrot_merger3, star4_vrot_merger3)
star1_AD_no_merger3 = va(gas_vrot_no_merger3, star1_vrot_no_merger3)
star2_AD_no_merger3 = va(gas_vrot_no_merger3, star2_vrot_no_merger3)
star3_AD_no_merger3 = va(gas_vrot_no_merger3, star3_vrot_no_merger3)
star4_AD_no_merger3 = va(gas_vrot_no_merger3, star4_vrot_no_merger3)

#get rid of nan values
star1_AD_merger1 = star1_AD_merger1[~np.isnan(star1_AD_merger1)]
star2_AD_merger1 = star2_AD_merger1[~np.isnan(star2_AD_merger1)]
star3_AD_merger1 = star3_AD_merger1[~np.isnan(star3_AD_merger1)]
star4_AD_merger1 = star4_AD_merger1[~np.isnan(star4_AD_merger1)]
star1_AD_no_merger1 = star1_AD_no_merger1[~np.isnan(star1_AD_no_merger1)]
star2_AD_no_merger1 = star2_AD_no_merger1[~np.isnan(star2_AD_no_merger1)]
star3_AD_no_merger1 = star3_AD_no_merger1[~np.isnan(star3_AD_no_merger1)]
star4_AD_no_merger1 = star4_AD_no_merger1[~np.isnan(star4_AD_no_merger1)]

star1_AD_merger2 = star1_AD_merger2[~np.isnan(star1_AD_merger2)]
star2_AD_merger2 = star2_AD_merger2[~np.isnan(star2_AD_merger2)]
star3_AD_merger2 = star3_AD_merger2[~np.isnan(star3_AD_merger2)]
star4_AD_merger2 = star4_AD_merger2[~np.isnan(star4_AD_merger2)]
star1_AD_no_merger2 = star1_AD_no_merger2[~np.isnan(star1_AD_no_merger2)]
star2_AD_no_merger2 = star2_AD_no_merger2[~np.isnan(star2_AD_no_merger2)]
star3_AD_no_merger2 = star3_AD_no_merger2[~np.isnan(star3_AD_no_merger2)]
star4_AD_no_merger2 = star4_AD_no_merger2[~np.isnan(star4_AD_no_merger2)]

star1_AD_merger3 = star1_AD_merger3[~np.isnan(star1_AD_merger3)]
star2_AD_merger3 = star2_AD_merger3[~np.isnan(star2_AD_merger3)]
star3_AD_merger3 = star3_AD_merger3[~np.isnan(star3_AD_merger3)]
star4_AD_merger3 = star4_AD_merger3[~np.isnan(star4_AD_merger3)]
star1_AD_no_merger3 = star1_AD_no_merger3[~np.isnan(star1_AD_no_merger3)]
star2_AD_no_merger3 = star2_AD_no_merger3[~np.isnan(star2_AD_no_merger3)]
star3_AD_no_merger3 = star3_AD_no_merger3[~np.isnan(star3_AD_no_merger3)]
star4_AD_no_merger3 = star4_AD_no_merger3[~np.isnan(star4_AD_no_merger3)]

#saving AD data from merger group 1 to use in comp_fig.py
star1_data_merger = calc_errors(star1_AD_merger1)
star2_data_merger = calc_errors(star2_AD_merger1)
star3_data_merger = calc_errors(star3_AD_merger1)
star4_data_merger = calc_errors(star4_AD_merger1)
star1_data_no_merger = calc_errors(star1_AD_no_merger1)
star2_data_no_merger = calc_errors(star2_AD_no_merger1)
star3_data_no_merger = calc_errors(star3_AD_no_merger1)
star4_data_no_merger = calc_errors(star4_AD_no_merger1)
median_mergers = np.array([star1_data_merger[0], star2_data_merger[0], star3_data_merger[0], star4_data_merger[0]])
lower_error_merger = np.array([star1_data_merger[1], star2_data_merger[1], star3_data_merger[1], star4_data_merger[1]])
upper_error_merger = np.array([star1_data_merger[2], star2_data_merger[2], star3_data_merger[2], star4_data_merger[2]])
median_no_mergers = np.array([star1_data_no_merger[0],  star2_data_no_merger[0], star3_data_no_merger[0], star4_data_no_merger[0]])
lower_error_no_merger = np.array([star1_data_no_merger[1], star2_data_no_merger[1], star3_data_no_merger[1], star4_data_no_merger[1]])
upper_error_no_merger = np.array([star1_data_no_merger[2], star2_data_no_merger[2], star3_data_no_merger[2], star4_data_no_merger[2]])
np.savetxt('/Volumes/Titan/analogs/TNGdata/comp_fig_data.txt', np.c_[median_mergers, lower_error_merger, upper_error_merger, median_no_mergers, lower_error_no_merger, upper_error_no_merger], header='median AD merger, lower error, upper error, median AD no mergers, lower error, upper error')
# #plots ===========================================================================================================================
# def single_plot():
# 	rc('font', family = 'serif')
# 	fig, ax=plt.subplots(1)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	plt.tick_params(which='both', width=1)
# 	plt.tick_params(which='major', length=7)
# 	plt.tick_params(which='minor', length=4)
# 	plt.tick_params(labelsize=12) 
# 	plt.minorticks_on()

# #AD histograms
# f, axes= plt.subplots(3,2, sharey=True, sharex=True, figsize=(16,18))

# axes[0,0].hist(star1_AD_merger1, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[0,0].hist(star2_AD_merger1, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[0,0].hist(star3_AD_merger1, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[0,0].hist(star4_AD_merger1, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[0,1].hist(star1_AD_no_merger1, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_no_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[0,1].hist(star2_AD_no_merger1, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_no_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[0,1].hist(star3_AD_no_merger1, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_no_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[0,1].hist(star4_AD_no_merger1, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_no_merger1),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[0,0].legend(loc=1, frameon=False, fontsize=10)
# axes[0,1].legend(loc=1, frameon=False)
# axes[0,0].annotate('Recent (4 Gyrs) 4:1 Merger', xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[0,0].annotate('{} Halos'.format(num_mergers1), xy=(-105, 0.054), horizontalalignment='left', fontsize=13)
# axes[0,1].annotate('No Recent 4:1 Merger', xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[0,1].annotate('{} Halos'.format(num_no_mergers1), xy=(-105, 0.054), horizontalalignment='left', fontsize=13)

# axes[1,0].hist(star1_AD_merger2, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[1,0].hist(star2_AD_merger2, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[1,0].hist(star3_AD_merger2, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[1,0].hist(star4_AD_merger2, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[1,1].hist(star1_AD_no_merger2, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_no_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[1,1].hist(star2_AD_no_merger2, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_no_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[1,1].hist(star3_AD_no_merger2, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_no_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[1,1].hist(star4_AD_no_merger2, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_no_merger2),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[1,0].legend(loc=1, frameon=False, fontsize=10)
# axes[1,1].legend(loc=1, frameon=False)
# axes[1,0].annotate('Recent ({} Gyrs) 4:1'.format(round(np.median(med_time),1)), xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[1,0].annotate('Merger', xy=(-99, 0.054), horizontalalignment='left', fontsize=13)
# axes[1,0].annotate('{} Halos'.format(num_mergers2), xy=(-105, 0.051), horizontalalignment='left', fontsize=13)
# axes[1,1].annotate('No Recent 4:1 Merger', xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[1,1].annotate('{} Halos'.format(num_no_mergers2), xy=(-105, 0.054), horizontalalignment='left', fontsize=13)

# axes[2,0].hist(star1_AD_merger3, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[2,0].hist(star2_AD_merger3, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[2,0].hist(star3_AD_merger3, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[2,0].hist(star4_AD_merger3, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[2,1].hist(star1_AD_no_merger3, bins=range(-100, 150, 10), label=r'$\rm \leq 1\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_AD_no_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# axes[2,1].hist(star2_AD_no_merger3, bins=range(-100, 150, 10), label=r'$\rm 1-5\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_AD_no_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# axes[2,1].hist(star3_AD_no_merger3, bins=range(-100, 150, 10), label=r'$\rm 5-10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_AD_no_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# axes[2,1].hist(star4_AD_no_merger3, bins=range(-100, 150, 10), label=r'$\rm \geq 10\ Gyr,\ \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_AD_no_merger3),2)) + r'$\rm \ km \ s^{-1}$', normed=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# axes[2,0].legend(loc=1, frameon=False, fontsize=10)
# axes[2,1].legend(loc=1, frameon=False)
# axes[2,0].annotate('Recent (12 Gyrs) 4:1 Merger', xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[2,0].annotate('{} Halos'.format(num_mergers3), xy=(-105, 0.054), horizontalalignment='left', fontsize=13)
# axes[2,1].annotate('No Recent 4:1 Merger', xy=(-105, 0.057), horizontalalignment='left', fontsize=13)
# axes[2,1].annotate('{} Halos'.format(num_no_mergers3), xy=(-105, 0.054), horizontalalignment='left', fontsize=13)

# for ax in axes[0,:]:
# 	ax.set_xlim(-115, 150)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	#ax.set_xticklabels([])

# for ax in axes[1,:]:
# 	ax.set_xlim(-115, 150)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	#ax.set_xticklabels([])

# for ax in axes[2,:]:
# 	ax.set_xlim(-115, 150)
# 	ax.set_xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=16)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	nbins = 6 
# 	ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# for ax in axes[:,0]:
# 	ax.set_ylabel(r'$ \rm PDF$', fontsize=16)
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	nbins = 8
# 	ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# for ax in axes[:,1]:
# 	#ax.set_ylabel(r'$ \rm PDF$', fontsize=16)
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# plt.subplots_adjust(wspace=.0, hspace=0)
# #plt.legend(frameon=False, fontsize=13)
# plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=16)
# plt.savefig('/Users/amandaquirk/Desktop/AD_merger_panels.png', bbox_inches='tight')
# plt.close()

# #RCs
# f, axes= plt.subplots(3,2, sharey=True, sharex=True, figsize=(16,18))

# axes[0,0].errorbar(r_bins, gas_vrot_merger1, yerr = gas_error_merger1, errorevery=20, c = 'darkgrey', label='gas')
# axes[0,0].errorbar(r_bins, star1_vrot_merger1, yerr = star1_error_merger1, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[0,0].errorbar(r_bins, star2_vrot_merger1, yerr = star2_error_merger1, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[0,0].errorbar(r_bins, star3_vrot_merger1, yerr = star3_error_merger1, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[0,0].errorbar(r_bins, star4_vrot_merger1, yerr = star4_error_merger1, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[0,1].errorbar(r_bins, gas_vrot_no_merger1,yerr = gas_error_no_merger1, errorevery=20, c = 'darkgrey', label='gas')
# axes[0,1].errorbar(r_bins, star1_vrot_no_merger1, yerr = star1_error_no_merger1, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[0,1].errorbar(r_bins, star2_vrot_no_merger1, yerr = star2_error_no_merger1, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[0,1].errorbar(r_bins, star3_vrot_no_merger1, yerr = star3_error_no_merger1, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[0,1].errorbar(r_bins, star4_vrot_no_merger1, yerr = star4_error_no_merger1, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[0,0].legend(loc=2, frameon=False, fontsize=13)
# axes[0,0].annotate('Recent (4 Gyrs) 4:1 Merger', xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[0,1].annotate('No Recent 4:1 Merger', xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[0,0].annotate('{} Halos'.format(num_mergers1), xy=(1, 230), horizontalalignment='left', fontsize=13)
# axes[0,1].annotate('{} Halos'.format(num_no_mergers1), xy=(1, 230), horizontalalignment='left', fontsize=13)

# axes[1,0].errorbar(r_bins, gas_vrot_merger2, yerr = gas_error_merger2, errorevery=20, c = 'darkgrey', label='gas')
# axes[1,0].errorbar(r_bins, star1_vrot_merger2, yerr = star1_error_merger2, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[1,0].errorbar(r_bins, star2_vrot_merger2, yerr = star2_error_merger2, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[1,0].errorbar(r_bins, star3_vrot_merger2, yerr = star3_error_merger2, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[1,0].errorbar(r_bins, star4_vrot_merger2, yerr = star4_error_merger2, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[1,1].errorbar(r_bins, gas_vrot_no_merger2,yerr = gas_error_no_merger2, errorevery=20, c = 'darkgrey', label='gas')
# axes[1,1].errorbar(r_bins, star1_vrot_no_merger2, yerr = star1_error_no_merger2, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[1,1].errorbar(r_bins, star2_vrot_no_merger2, yerr = star2_error_no_merger2, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[1,1].errorbar(r_bins, star3_vrot_no_merger2, yerr = star3_error_no_merger2, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[1,1].errorbar(r_bins, star4_vrot_no_merger2, yerr = star4_error_no_merger2, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[1,0].legend(loc=2, frameon=False, fontsize=13)
# axes[1,0].annotate('Recent ({} Gyrs) 4:1'.format(round(np.median(med_time),1)), xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[1,1].annotate('No Recent 4:1 Merger', xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[1,0].annotate('{} Halos'.format(num_mergers2), xy=(1, 230), horizontalalignment='left', fontsize=13)
# axes[1,1].annotate('{} Halos'.format(num_no_mergers2), xy=(1, 230), horizontalalignment='left', fontsize=13)

# axes[2,0].errorbar(r_bins, gas_vrot_merger3, yerr = gas_error_merger3, errorevery=20, c = 'darkgrey', label='gas')
# axes[2,0].errorbar(r_bins, star1_vrot_merger3, yerr = star1_error_merger3, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[2,0].errorbar(r_bins, star2_vrot_merger3, yerr = star2_error_merger3, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[2,0].errorbar(r_bins, star3_vrot_merger3, yerr = star3_error_merger3, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[2,0].errorbar(r_bins, star4_vrot_merger3, yerr = star4_error_merger3, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[2,1].errorbar(r_bins, gas_vrot_no_merger3,yerr = gas_error_no_merger3, errorevery=20, c = 'darkgrey', label='gas')
# axes[2,1].errorbar(r_bins, star1_vrot_no_merger3, yerr = star1_error_no_merger3, errorevery=21, c = 'b', linestyle='--', alpha = 0.6, label=r'$\rm \leq 1\ Gyr$')
# axes[2,1].errorbar(r_bins, star2_vrot_no_merger3, yerr = star2_error_no_merger3, errorevery=19, c = 'm', alpha = 0.6, label=r'$\rm 1-5\ Gyr$')
# axes[2,1].errorbar(r_bins, star3_vrot_no_merger3, yerr = star3_error_no_merger3, errorevery=18, c = 'green', alpha = 0.6, label=r'$\rm 5-10\ Gyr$')
# axes[2,1].errorbar(r_bins, star4_vrot_no_merger3, yerr = star4_error_no_merger3, errorevery=17, c = 'r', linestyle=':', alpha = 0.6, label=r'$\rm \geq 10\ Gyr$')
# axes[2,0].legend(loc=2, frameon=False, fontsize=13)
# axes[2,0].annotate('Recent (12 Gyrs) 4:1 Merger', xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[2,1].annotate('No Recent 4:1 Merger', xy=(1, 240), horizontalalignment='left', fontsize=13)
# axes[2,0].annotate('{} Halos'.format(num_mergers3), xy=(1, 230), horizontalalignment='left', fontsize=13)
# axes[2,1].annotate('{} Halos'.format(num_no_mergers3), xy=(1, 230), horizontalalignment='left', fontsize=13)

# for ax in axes[0,:]:
# 	ax.set_xlim(0, 20)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	#ax.set_xticklabels([])

# for ax in axes[1,:]:
# 	ax.set_xlim(0, 20)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# 	#ax.set_xticklabels([])

# for ax in axes[2,:]:
# 	ax.set_xlim(0, 20)
# 	ax.set_xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=16)
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	ax.legend(fontsize=11.5, frameon=False)
# 	nbins = 6 
# 	ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# for ax in axes[:,0]:
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	nbins = 8
# 	ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
# 	ax.minorticks_on()
# 	ax.set_ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=16)
# 	ax.set_ylim(0, 270)
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# for ax in axes[:,1]:
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=12) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)

# plt.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/rc_merger_panels.png', bbox_inches='tight')
# plt.close()









