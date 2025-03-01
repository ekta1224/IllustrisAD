import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc

median_merger, lerror_merger, uerror_merger, median_no, lerror_no, uerror_no = np.loadtxt('/Volumes/Titan/analogs/TNGdata/comp_fig_data.txt', unpack=True)

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

m31_age = [30e6, 400e6, 2e9, 4e9]
m31_ad = [-8.15, 17.69, 50.43, 62.97]
m31_upper_err = [0.74, 3.29, 1.09, 0.59]
m31_lower_err = [0.72, 2.78, 1.26, 0.40]

m31_x_err = np.array(m31_age) / 2
#print(m31_x_upper_err)
#print(m31_x_lower_err)

illustris_age = [0, 1e9, 1e9, 5e9, 5e9, 10e9, 10e9, 13.7e9] 

merger_ad =        [median_merger[0], median_merger[0], median_merger[1], median_merger[1], median_merger[2], median_merger[2], median_merger[3], median_merger[3]]
merger_upper_err = [lerror_merger[0], lerror_merger[0], lerror_merger[1], lerror_merger[1], lerror_merger[2], lerror_merger[2], lerror_merger[3], lerror_merger[3]]
merger_lower_err = [uerror_merger[0], uerror_merger[0], uerror_merger[1], uerror_merger[1], uerror_merger[2], uerror_merger[2], uerror_merger[3], uerror_merger[3]] 

no_merger_ad =        [median_no[0], median_no[0], median_no[1], median_no[1], median_no[2], median_no[2], median_no[3], median_no[3]]
no_merger_upper_err = [lerror_no[0], lerror_no[0], lerror_no[1], lerror_no[1], lerror_no[2], lerror_no[2], lerror_no[3], lerror_no[3]]
no_merger_lower_err = [uerror_no[0], uerror_no[0], uerror_no[1], uerror_no[1], uerror_no[2], uerror_no[2], uerror_no[3], uerror_no[3]]

single_plot()
plt.plot([illustris_age[0], illustris_age[1]], [merger_ad[0],merger_ad[1]], color='cornflowerblue', linestyle='--')
plt.plot([illustris_age[2], illustris_age[3]], [merger_ad[2],merger_ad[3]], color='cornflowerblue', linestyle='--')
plt.plot([illustris_age[4], illustris_age[5]], [merger_ad[4],merger_ad[5]], color='cornflowerblue', linestyle='--')
plt.plot([illustris_age[6], illustris_age[7]], [merger_ad[6],merger_ad[7]], color='cornflowerblue', linestyle='--')
plt.fill_between(illustris_age, np.array(merger_ad) + np.array(merger_upper_err), np.array(merger_ad) - np.array(merger_lower_err), color='cornflowerblue', alpha=0.5)
plt.plot([illustris_age[0], illustris_age[1]], [no_merger_ad[0], no_merger_ad[1]], color='hotpink', linestyle='--')
plt.plot([illustris_age[2], illustris_age[3]], [no_merger_ad[2], no_merger_ad[3]], color='hotpink', linestyle='--')
plt.plot([illustris_age[4], illustris_age[5]], [no_merger_ad[4], no_merger_ad[5]], color='hotpink', linestyle='--')
plt.plot([illustris_age[6], illustris_age[7]], [no_merger_ad[6], no_merger_ad[7]], color='hotpink', linestyle='--')
plt.fill_between(illustris_age, np.array(no_merger_ad) + np.array(no_merger_upper_err), np.array(no_merger_ad) - np.array(no_merger_lower_err), color='hotpink', alpha=0.5)
plt.plot([4e9, 4e9], [-40, 140], c='grey', linestyle='-')
eb = plt.errorbar(m31_age, m31_ad, yerr=[m31_lower_err, m31_upper_err], xerr=[m31_x_err, m31_x_err], fmt='o', c='k', label=r'$\rm M31\ Observations$')
eb[-1][0].set_linestyle(':')
plt.xlabel(r'$ \rm Stellar\ Age\ (yr)$', fontsize=16)
plt.ylabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=16)
plt.ylim(-35, 150)
plt.xscale('log')
plt.legend(loc=2, frameon=False, fontsize=13)
plt.annotate(r'$\rm Quirk\ et\ al.\ 2019$', (5.55e7, 108), fontsize=13)
plt.annotate(r'$\rm Analogs\ with\ a\ 4:1\ Merger$', (3e7, 95), fontsize=13, color='cornflowerblue')
plt.annotate(r'$\rm in\ the\ last\ 4\ Gyr$', (3e7, 85), fontsize=13, color='cornflowerblue')
plt.annotate(r'$\rm Analogs\ without\ a\ 4:1\ Merger$', (3e7, 76), fontsize=13, color='hotpink')
plt.annotate(r'$\rm in\ the\ last\ 4\ Gyr$', (3e7, 66), fontsize=13, color='hotpink')
#plt.show()
plt.savefig('/Users/amandaquirk/Desktop/M31_comp.png', bbox_inches='tight')



