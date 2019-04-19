import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 

'''
looks at median AD vs halo properties: first creates a data file, then does some plotting
'''

# save data for each analog to a file ==================================================================
#read in per
halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt') #halo IDs
ID, Mmerger_times = np.loadtxt('../data/M31analogs_halo_props.txt', usecols=(0,9,), unpack=True) #halo properties file

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

group1_AD = np.zeros_like(halos)
group2_AD = np.zeros_like(halos)
group3_AD = np.zeros_like(halos)
group4_AD = np.zeros_like(halos)
group1_error = np.zeros_like(halos)
group2_error = np.zeros_like(halos)
group3_error = np.zeros_like(halos)
group4_error = np.zeros_like(halos)
time_M_merger = np.zeros_like(halos)
for i in range(len(halos)):
	#print(halos[i])

	#the halo proos file has all of the analogs Ekta found, so need to match IDs to get the merger info
	N = np.where(halos[i] == ID)

	#rotation velocity indo
	gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/FRIEND/analogs/data/{}_vrot_gas_smooted.txt'.format(int(halos[i])), usecols=(1,2,3,4,5,), unpack=True) #km/s

	#calculate the median AD for each group across all of the radial bins
	group1_AD[i] = np.nanmedian(va(gas_vrot, star1_vrot))
	group2_AD[i] = np.nanmedian(va(gas_vrot, star2_vrot))
	group3_AD[i] = np.nanmedian(va(gas_vrot, star3_vrot))
	group4_AD[i] = np.nanmedian(va(gas_vrot, star4_vrot))
	group1_error[i] = np.nanstd(va(gas_vrot, star1_vrot))
	#some of group1 didn't have enough data, so using below two lines to remove the halos from the whole anaylsis
	if group1_error[i] == 0:
		print(halos[i])
	group2_error[i] = np.nanstd(va(gas_vrot, star2_vrot))
	group3_error[i] = np.nanstd(va(gas_vrot, star3_vrot))
	group4_error[i] = np.nanstd(va(gas_vrot, star4_vrot))
	time_M_merger[i] = Mmerger_times[N]

np.savetxt('/Volumes/FRIEND/analogs/data/AD_halo_props.txt', np.c_[halos, group1_AD, group2_AD, group3_AD, group4_AD, group1_error, group2_error, group3_error, group4_error, time_M_merger], fmt="%-s", delimiter='\t', header='ID, star 1 group AD, star 2 group AD, star 3 group AD, star 4 group AD, star 1 group AD error, star 2 group AD error, star 3 group AD error, star 4 group AD error, time since last 4:1 merger') 
# ======================================================================================================
#plotting 
def single_plot():
	rc('font', family = 'serif')
	fig, ax=plt.subplots(1)
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()

#median AD
#doing linear fits 
times = np.linspace(0, 13.8, 50)
m1, b1 = np.polyfit(time_M_merger, group1_AD, 1, w=(1/group1_error))
m2, b2 = np.polyfit(time_M_merger, group2_AD, 1, w=(1/group2_error))
m3, b3 = np.polyfit(time_M_merger, group3_AD, 1, w=(1/group3_error))
m4, b4 = np.polyfit(time_M_merger, group4_AD, 1, w=(1/group4_error))

single_plot()
plt.scatter(time_M_merger, group1_AD, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1, y={}x+{}'. format(round(m1,3), round(b1,2)))
plt.scatter(time_M_merger, group2_AD, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2, y={}x+{}'. format(round(m2,3), round(b2,2)))
plt.scatter(time_M_merger, group3_AD, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3, y={}x+{}'. format(round(m3,3), round(b3,2)))
plt.scatter(time_M_merger, group4_AD, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4, y={}x+{}'. format(round(m4,3), round(b4,2)))
plt.plot(times, m1 * times + b1, c='b')
plt.plot(times, m2 * times + b2, c='m')
plt.plot(times, m3 * times + b3, c='green')
plt.plot(times, m4 * times + b4, c='r')
plt.xlabel('Time Since Last 4:1 Merger (Gyr)')
plt.ylabel('Median AD (km/s)')
plt.ylim(-50,)
plt.legend(frameon= False)
plt.savefig('/Volumes/FRIEND/analogs/plots/merger/AD_major_merger.png', bbox_inches='tight')
plt.close()


