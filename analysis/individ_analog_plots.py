import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc

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

analog = 419510

r, gas_vrot, star1_vrot, star2_vrot, star3_vrot, star4_vrot = np.loadtxt('/Volumes/Titan/analogs/data/{}_vrot_gas_smooted.txt'.format(int(analog)), unpack=True)

def va(v_gas, v_star):
		return v_gas - v_star
	
star1_ad = va(gas_vrot, star1_vrot)
star2_ad = va(gas_vrot, star2_vrot)
star3_ad = va(gas_vrot, star3_vrot)
star4_ad = va(gas_vrot, star4_vrot)

#remove nan values (that resulted from having empty radial bins)
star1_ad = star1_ad[~np.isnan(star1_ad)]
star2_ad = star2_ad[~np.isnan(star2_ad)]
star3_ad = star3_ad[~np.isnan(star3_ad)]
star4_ad = star4_ad[~np.isnan(star4_ad)]

#rc
# single_plot()
# plt.scatter(r, gas_vrot, c = 'darkgrey', s=12, label='gas')
# plt.scatter(r, star1_vrot, c = 'b', alpha = 0.6, s=10, marker='^', label=r'$ \rm \leq 1\ Gyr$')
# plt.scatter(r, star2_vrot, c = 'm', alpha = 0.6, s=10, marker='P', label=r'$ \rm 1-5\ Gyr$')
# plt.scatter(r, star3_vrot, c = 'green', alpha = 0.6, s=10, marker='s', label=r'$ \rm 5-10\ Gyr$')
# plt.scatter(r, star4_vrot, c = 'r', alpha = 0.6, s=14, marker='_', label=r'$ \rm \geq 10\ Gyr$')
# plt.ylim(0,250)
# plt.xlim(2,20)
# plt.ylabel(r'$ \rm Rotation\ Velocity:\ \itv_{\rm rot}\ \rm(km\ s^{-1})$', fontsize=13)
# plt.xlabel(r'$\rm Radial\ Distance:\ \it r\ \rm (kpc)$', fontsize=13)
# plt.legend(loc=2, frameon=False)
# plt.savefig('/Volumes/Titan/analogs/plots/paperplots/{}_rc.pdf'.format(int(analog)), bbox_inches='tight') 
# plt.close()

#ad hist
single_plot()
plt.hist(star1_ad, bins=range(-120, 160, 10), label=r'$\rm \leq 1\ Gyr, \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star1_ad),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(star2_ad, bins=range(-120, 160, 10), label=r'$\rm 1-5\ Gyr, \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star2_ad),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
plt.hist(star3_ad, bins=range(-120, 160, 10), label=r'$\rm 5-10\ Gyr, \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star3_ad),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
plt.hist(star4_ad, bins=range(-120, 160, 10), label=r'$\rm \geq 10\ Gyr, \overline{v_{a}}$' + r'$={}$'.format(round(np.median(star4_ad),2)) + r'$\rm \ km \ s^{-1}$', density=True, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
plt.legend(loc=1, frameon=False)
plt.xlim(-120,180)
plt.ylabel(r'$ \rm PDF$', fontsize=16)
plt.xlabel(r'$ \rm Asymmetric\ Drift:\ \itv_{a}\ \rm(km\ s^{-1})$', fontsize=13)
plt.savefig('/Users/amandaquirk/Desktop/{}_ad.pdf'.format(int(analog)), bbox_inches='tight')

#spatial plot