import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc

halos = np.loadtxt('../data/M31analog_IDs_IllustrisAD.txt')
ID, mass_tot, m_star, Merger_times, num_M, num_m, ex_situ_mass, in_situ_mass, accret_mass = np.loadtxt('../data/M31analogs_halo_props.txt', usecols=(0, 7, 9, 10, 11, 12, 13, 14, 15), unpack=True) #halo properties file
h = 0.704

inds = []
for i in range(len(halos)):
	#print(int(halos[i]))

	#the halo proos file has all of the analogs Ekta found, so need to match IDs to get the merger info
	N = np.where(halos[i] == ID)
	inds.append(N[0][0])

Mvir = mass_tot[inds] * 10**10 / h #Msun
Mstar = m_star[inds] * 10**10 / h #Msun
M_time = Merger_times[inds]
num_Major = num_M[inds]
num_minor = num_m[inds]
num_tot = num_Major + num_minor
id_prim = ID[inds]

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
	plt.tick_params(labelsize=13) 
	plt.minorticks_on()

# #merger time histogram
# single_plot()
# plt.hist(Merger_times, color='k', fill=False, histtype='step', linewidth=2, label=r'$\rm 249\ Analogs\ $')
# plt.hist(M_time, color='k', fill=False,  hatch='//', histtype='step', linewidth=2, label=r'$\rm Primary\ Sample\ of\ 150\ Analogs\ $')
# plt.plot([9.4, 9.4], [0,60], color='darkgrey', linestyle='--', linewidth=3)
# plt.plot([np.median(M_time), np.median(M_time)], [0,60], color='darkgrey', linestyle='--', linewidth=3)
# plt.ylim(0,55)
# plt.xlabel(r'$\rm Time\ Since\ Last\ 4:1\ Merger\ (Gyr)$', fontsize=16)
# plt.ylabel(r'$\rm N$', fontsize=16)
# plt.legend(frameon=False, fontsize=11.5, loc=2)
# plt.savefig('/Users/amandaquirk/Desktop/merger_times.png')
# plt.close()

#mvir and mstar scatter plot
# single_plot()
# plt.scatter(Mvir, Mstar, c='k')
# plt.xlabel(r'$M_{vir} \ (M_{sun})$', fontsize=13)
# plt.ylabel(r'$M_{\star}\ (M_{sun})$', fontsize=13)
# plt.savefig('/Users/amandaquirk/Desktop/halo_masses_scatter.png')
# plt.close()

# mvir and mstar histogram plot
# f, axes= plt.subplots(1,2, sharey=True, sharex=False, figsize=(15,6))
# for ax in axes:
# 	ax.tick_params(axis='x',which='both',bottom='on',top='off', direction='out')
# 	ax.tick_params(axis='x',which='both',top='on', direction='in')
# 	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
# 	ax.tick_params(axis='y',which='both',right='on', direction='in')
# 	ax.tick_params(which='both', width=1)
# 	ax.tick_params(which='major', length=7)
# 	ax.tick_params(which='minor', length=4)
# 	ax.tick_params(labelsize=13) 
# 	ax.minorticks_on()
# 	for axis in ['top','bottom','left','right']:
# 	        ax.spines[axis].set_linewidth(1)
# axes[0].hist(mass_tot * 10**10 / h , color='k', fill=False, histtype='step', linewidth=2, label=r'$\rm 249\ Analogs\ $')	        
# axes[0].hist(Mvir, color='k', fill=False,  hatch='//', histtype='step', linewidth=2, label=r'$\rm Primary\ Sample\ of\ 150\ Analogs\ $')
# axes[1].hist(m_star * 10**10 / h , color='k', fill=False, histtype='step', linewidth=2, label=r'$\rm 249\ Analogs\ $')
# axes[1].hist(Mstar, color='k', fill=False,  hatch='//', histtype='step', linewidth=2, label=r'$\rm Primary\ Sample\ of\ 150\ Analogs\ $')
# axes[0].set_xlabel(r'$M_{vir} \ (M_{sun})$', fontsize=16)
# axes[1].set_xlabel(r'$M_{\star} \ (M_{sun})$', fontsize=16)
# axes[0].set_ylabel(r'$\rm N$', fontsize=16)
# axes[0].legend(frameon=False, fontsize=11.5)
# axes[1].legend(frameon=False, fontsize=11.5)
# plt.subplots_adjust(wspace=.0, hspace=0)
# plt.savefig('/Users/amandaquirk/Desktop/halo_masses_hist.png')
# plt.close()

#number of merger histogram
# single_plot()
# plt.hist(num_Major, bins=range(0,14, 1), fill=False, histtype='step', linewidth=2, hatch='//', label=r'$\rm 4:1\ Mergers$')
# plt.hist(num_minor, bins=range(0,14, 1), color='r', alpha=.6, fill=False, histtype='step', linewidth=2, label=r'$\rm 10:1\ Mergers$')
# plt.hist(num_tot, bins=range(0,14, 1), color='k', fill=False, histtype='step', linewidth=2, label=r'$\rm Total\ Mergers$')
# plt.xlabel(r'$\rm Number\ of\ Mergers$', fontsize=16)
# plt.ylabel(r'$\rm N$', fontsize=16)
# plt.legend(frameon=False, fontsize=11.5)
# plt.savefig('/Users/amandaquirk/Desktop/num_mergers.png')
# plt.close()








