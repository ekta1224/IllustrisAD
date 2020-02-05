import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 

'''
looks at median AD vs halo properties: first creates a data file, then does some plotting
'''

# save data for each analog to a file ==================================================================
#read in per
halos = np.loadtxt('../data/M31analog_noM33_IDs.txt') #halo IDs
ID, mass_tot, Mmerger_times, num_M, num_m, ex_situ_mass, in_situ_mass, accret_mass = np.loadtxt('../data/M31analogs_halo_props_strictly_noM33.txt', usecols=(0, 9, 10, 11, 12, 13, 14, 15), unpack=True) #halo properties file

#Illustris value
h = 0.704

#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

#below was used to get the median AD values for each analog. no need to rereun unless change to Illustris_curves.py
group1_AD = np.zeros_like(halos)
group2_AD = np.zeros_like(halos)
group3_AD = np.zeros_like(halos)
group4_AD = np.zeros_like(halos)
group1_error = np.zeros_like(halos)
group2_error = np.zeros_like(halos)
group3_error = np.zeros_like(halos)
group4_error = np.zeros_like(halos)
time_M_merger = np.zeros_like(halos)
num_M_merger = np.zeros_like(halos)
num_m_merger = np.zeros_like(halos)
ex_mass = np.zeros_like(halos)
in_mass = np.zeros_like(halos)
accreted_mass = np.zeros_like(halos)
stellar_mass = np.zeros_like(halos)
for i in range(len(halos)):
	#print(int(halos[i]))

	#the halo proos file has all of the analogs Ekta found, so need to match IDs to get the merger info
	N = np.where(halos[i] == ID)
	#print(N)

	#rotation velocity indo
	gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/FRIEND/analogs/data/{}_vrot_gas_smoothed.txt'.format(int(halos[i])), usecols=(1,2,3,4,5,), unpack=True) #km/s

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
	num_M_merger[i] = num_M[N]
	num_m_merger[i] = num_m[N]
	stellar_mass[i] = mass_tot[N] * 10**10 / h #Msun
	ex_mass[i] = ex_situ_mass[N] * 10**10 / h #Msun
	in_mass[i] = in_situ_mass[N] * 10**10 / h #Msun
	accreted_mass[i] = accret_mass[N] * 10**10 / h #Msun

np.savetxt('/Volumes/FRIEND/analogs/data/AD_halo_props_noM33.txt', np.c_[halos, group1_AD, group2_AD, group3_AD, group4_AD, group1_error, group2_error, group3_error, group4_error, stellar_mass, time_M_merger, num_M_merger, num_m_merger, ex_mass, in_mass, accreted_mass], fmt="%-s", delimiter='\t', header='ID, star 1 group AD, star 2 group AD, star 3 group AD, star 4 group AD, star 1 group AD error, star 2 group AD error, star 3 group AD error, star 4 group AD error, total stellar mass (Msun), time since last 4:1 merger, number of 4:1 mergers, number of 10:1 mergers, mass formed ex situ (Msun), mass formed in situ (Msun), mass accreted (Msun)') 
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

group1_AD, group2_AD, group3_AD, group4_AD, group1_error, group2_error, group3_error, group4_error, stellar_mass, time_M_merger, num_M_merger, num_m_merger, ex_mass, in_mass, accreted_mass = np.loadtxt('/Volumes/FRIEND/analogs/data/AD_halo_props_noM33.txt', usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), unpack=True)

def med_AD_plot(parameter, label):

	if label == 'time_M_merger':
		x_label = 'Time Since Last 4:1 Meger (Gyr)'
		plotname = 'AD_time_major_merger'
	if label =='num_M_merger':
		x_label = r'$\rm Number\ of\ 4:1\ Mergers$'
		plotname = 'AD_num_major_mergers'
	if label == 'num_m_merger':
		x_label = 'Number of 10:1 Mergers'
		plotname = 'AD_num_minor_mergers'
	if label == 'ex_mass_frac':
		x_label = r'$\rm Stellar\ Mass\ Fraction\ Formed\ Ex\ Situ$'
		plotname = 'AD_ex_situ_mass'
	if label == 'in_mass':
		x_label = 'Mass Formed In Situ (Msun)'
		plotname = 'AD_in_situ_mass'
	if label == 'accreted_mass':
		x_label = 'Mass Accreted from Completed Mergers (Msun)'
		plotname = 'AD_accreted_mass'

	#linear plots
	x_vals = np.linspace(min(parameter), max(parameter), 50)
	m1, b1 = np.polyfit(parameter, group1_AD, 1, w=(1/group1_error))
	m2, b2 = np.polyfit(parameter, group2_AD, 1, w=(1/group2_error))
	m3, b3 = np.polyfit(parameter, group3_AD, 1, w=(1/group3_error))
	m4, b4 = np.polyfit(parameter, group4_AD, 1, w=(1/group4_error))

	single_plot()
	plt.scatter(parameter, group1_AD, c = 'b', alpha = 0.8, s=25, marker='^', label=r'$\rm \leq 1\ Gyr$')#(round(m1,3), round(b1,2)))
	plt.scatter(parameter, group2_AD, c = 'm', alpha = 0.8, s=25, marker='P', label=r'$\rm 1-5\ Gyr$')#round(m2,3), round(b2,2)))
	plt.scatter(parameter, group3_AD, c = 'green', alpha = 0.8, s=25, marker='s', label=r'$\rm 5-10\ Gyr$')#(round(m3,3), round(b3,2)))
	plt.scatter(parameter, group4_AD, c = 'r', s=35, marker='_', label=r'$\rm \geq 10\ Gyr$')#round(m4,3), round(b4,2)))
	# plt.plot(x_vals, m1 * x_vals + b1, c='b')
	# plt.plot(x_vals, m2 * x_vals + b2, c='m')
	# plt.plot(x_vals, m3 * x_vals + b3, c='green')
	# plt.plot(x_vals, m4 * x_vals + b4, c='r')
	plt.xlabel('{}'.format(x_label), fontsize=13)
	plt.ylabel(r'$\rm Median\ AD\ (km\ s^{-1})$', fontsize=13)
	plt.ylim(-50,)
	#plt.legend(frameon= True, fontsize=10)
	plt.savefig('/Users/amandaquirk/Desktop/{}.pdf'.format(plotname), bbox_inches='tight')
	plt.close()

ex_mass_frac = ex_mass / stellar_mass
med_AD_plot(time_M_merger, 'time_M_merger')
med_AD_plot(num_M_merger, 'num_M_merger')
med_AD_plot(num_m_merger, 'num_m_merger')
med_AD_plot(ex_mass_frac, 'ex_mass_frac')
med_AD_plot(in_mass, 'in_mass')
med_AD_plot(accreted_mass, 'accreted_mass')


