import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 

group1_AD, group2_AD, group3_AD, group4_AD, group1_error, group2_error, group3_error, group4_error = np.loadtxt('/Volumes/FRIEND/analogs/data/AD_halo_props.txt', usecols=(1, 2, 3, 4, 5, 6, 7, 8,), unpack=True)

star1_med_disp_R, star2_med_disp_R, star3_med_disp_R, star4_med_disp_R, star1_med_disp_Z, star2_med_disp_Z, star3_med_disp_Z, star4_med_disp_Z, star1_med_disp_phi, star2_med_disp_phi, star3_med_disp_phi, star4_med_disp_phi = np.loadtxt('/Volumes/FRIEND/analogs/data/star_particle_dispersion.txt', usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), unpack=True)

z_r_1 = star1_med_disp_Z / star1_med_disp_R
phi_r_1 = star1_med_disp_phi / star1_med_disp_R

z_r_2 = star2_med_disp_Z / star2_med_disp_R
phi_r_2 = star2_med_disp_phi / star2_med_disp_R

z_r_3 = star3_med_disp_Z / star3_med_disp_R
phi_r_3 = star3_med_disp_phi / star3_med_disp_R

z_r_4 = star4_med_disp_Z / star4_med_disp_R
phi_r_4 = star4_med_disp_phi / star4_med_disp_R

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

# x_vals = np.linspace(-50, 200, 50)
# #plot one
# m1, b1 = np.polyfit(group1_AD, (star1_med_disp_phi**2), 1, w=(1/group1_error))
# m2, b2 = np.polyfit(group2_AD, (star2_med_disp_phi**2), 1, w=(1/group2_error))
# m3, b3 = np.polyfit(group3_AD, (star3_med_disp_phi**2), 1, w=(1/group3_error))
# m4, b4 = np.polyfit(group4_AD, (star4_med_disp_phi**2), 1, w=(1/group4_error))

# single_plot()
# plt.scatter(group1_AD, star1_med_disp_phi**2, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1, y={}x+{}'. format(m1, b1))#(round(m1,3), round(b1,2)))
# plt.scatter(group2_AD, star2_med_disp_phi**2, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2, y={}x+{}'. format(m2, b2))#round(m2,3), round(b2,2)))
# plt.scatter(group3_AD, star3_med_disp_phi**2, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3, y={}x+{}'. format(m3, b3))#(round(m3,3), round(b3,2)))
# plt.scatter(group4_AD, star4_med_disp_phi**2, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4, y={}x+{}'. format(m4, b4))#round(m4,3), round(b4,2)))
# plt.plot(x_vals, m1 * x_vals + b1, c='b')
# plt.plot(x_vals, m2 * x_vals + b2, c='m')
# plt.plot(x_vals, m3 * x_vals + b3, c='green')
# plt.plot(x_vals, m4 * x_vals + b4, c='r')
# plt.xlabel('Median AD (km/s)')
# plt.ylabel(r'$\sigma_{\phi}^{2}$')
# #plt.ylim(-50,)
# plt.legend(frameon= False)
# plt.savefig('/Users/amandaquirk/Desktop/med_AD_dsip_phi.png', bbox_inches='tight')
# plt.close()

# #plot two
# m1, b1 = np.polyfit(group1_AD, z_r_1, 1, w=(1/group1_error))
# m2, b2 = np.polyfit(group2_AD, z_r_2, 1, w=(1/group2_error))
# m3, b3 = np.polyfit(group3_AD, z_r_3, 1, w=(1/group3_error))
# m4, b4 = np.polyfit(group4_AD, z_r_4, 1, w=(1/group4_error))

# single_plot()
# plt.scatter(group1_AD, z_r_1, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1, y={}x+{}'. format(m1, b1))#(round(m1,3), round(b1,2)))
# plt.scatter(group2_AD, z_r_2, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2, y={}x+{}'. format(m2, b2))#round(m2,3), round(b2,2)))
# plt.scatter(group3_AD, z_r_3, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3, y={}x+{}'. format(m3, b3))#(round(m3,3), round(b3,2)))
# plt.scatter(group4_AD, z_r_4, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4, y={}x+{}'. format(m4, b4))#round(m4,3), round(b4,2)))
# plt.plot(x_vals, m1 * x_vals + b1, c='b')
# plt.plot(x_vals, m2 * x_vals + b2, c='m')
# plt.plot(x_vals, m3 * x_vals + b3, c='green')
# plt.plot(x_vals, m4 * x_vals + b4, c='r')
# plt.xlabel('Median AD (km/s)')
# plt.ylabel(r'$\sigma_{z} / \sigma_{R}$')
# #plt.ylim(-50,)
# plt.legend(frameon= False)
# plt.savefig('/Users/amandaquirk/Desktop/med_AD_dsip_z_r.png', bbox_inches='tight')
# plt.close()

# #plot 3
# m1, b1 = np.polyfit(group1_AD, phi_r_1, 1, w=(1/group1_error))
# m2, b2 = np.polyfit(group2_AD, phi_r_2, 1, w=(1/group2_error))
# m3, b3 = np.polyfit(group3_AD, phi_r_3, 1, w=(1/group3_error))
# m4, b4 = np.polyfit(group4_AD, phi_r_4, 1, w=(1/group4_error))

# single_plot()
# plt.scatter(group1_AD, phi_r_1, c = 'b', alpha = 0.6, s=10, marker='^', label='Group 1, y={}x+{}'. format(m1, b1))#(round(m1,3), round(b1,2)))
# plt.scatter(group2_AD, phi_r_2, c = 'm', alpha = 0.6, s=10, marker='P', label='Group 2, y={}x+{}'. format(m2, b2))#round(m2,3), round(b2,2)))
# plt.scatter(group3_AD, phi_r_3, c = 'green', alpha = 0.6, s=10, marker='s', label='Group 3, y={}x+{}'. format(m3, b3))#(round(m3,3), round(b3,2)))
# plt.scatter(group4_AD, phi_r_4, c = 'r', alpha = 0.6, s=14, marker='_', label='Group 4, y={}x+{}'. format(m4, b4))#round(m4,3), round(b4,2)))
# plt.plot(x_vals, m1 * x_vals + b1, c='b')
# plt.plot(x_vals, m2 * x_vals + b2, c='m')
# plt.plot(x_vals, m3 * x_vals + b3, c='green')
# plt.plot(x_vals, m4 * x_vals + b4, c='r')
# plt.xlabel('Median AD (km/s)')
# plt.ylabel(r'$\sigma_{\phi} / \sigma_{R}$')
# #plt.ylim(-50,)
# plt.legend(frameon= False)
# plt.savefig('/Users/amandaquirk/Desktop/med_AD_dsip_phi_r.png', bbox_inches='tight')
# plt.close()

# #histograms of velocity ellipsoid components
# single_plot()
# plt.hist(star1_med_disp_Z , bins=np.linspace(0, 150, 25), label='Group 1: med = {}'.format(np.median(star1_med_disp_Z)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# plt.hist(star2_med_disp_Z , bins=np.linspace(0, 150, 25), label='Group 2: med = {}'.format(np.median(star2_med_disp_Z)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# plt.hist(star3_med_disp_Z , bins=np.linspace(0, 150, 25), label='Group 3: med = {}'.format(np.median(star3_med_disp_Z)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# plt.hist(star4_med_disp_Z , bins=np.linspace(0, 150, 25), label='Group 4: med = {}'.format(np.median(star4_med_disp_Z)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# plt.legend(loc=2, frameon=False)
# plt.xlabel(r'$\sigma_{z}$')
# plt.savefig('/Users/amandaquirk/Desktop/med_z_hist.png', bbox_inches='tight')
# plt.close()

# single_plot()
# plt.hist(star1_med_disp_phi, bins=np.linspace(0, 150, 25), label='Group 1: med = {}'.format(np.median(star1_med_disp_phi)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# plt.hist(star2_med_disp_phi, bins=np.linspace(0, 150, 25), label='Group 2: med = {}'.format(np.median(star2_med_disp_phi)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# plt.hist(star3_med_disp_phi, bins=np.linspace(0, 150, 25), label='Group 3: med = {}'.format(np.median(star3_med_disp_phi)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# plt.hist(star4_med_disp_phi, bins=np.linspace(0, 150, 25), label='Group 4: med = {}'.format(np.median(star4_med_disp_phi)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# plt.legend(loc=2, frameon=False)
# plt.xlabel(r'$\sigma_{\phi}$')
# plt.savefig('/Users/amandaquirk/Desktop/med_phi_hist.png', bbox_inches='tight')
# plt.close()

# single_plot()
# plt.hist(star1_med_disp_R, bins=np.linspace(0, 150, 25), label='Group 1: med = {}'.format(np.median(star1_med_disp_R)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
# plt.hist(star2_med_disp_R, bins=np.linspace(0, 150, 25), label='Group 2: med = {}'.format(np.median(star2_med_disp_R)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
# plt.hist(star3_med_disp_R, bins=np.linspace(0, 150, 25), label='Group 3: med = {}'.format(np.median(star3_med_disp_R)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
# plt.hist(star4_med_disp_R, bins=np.linspace(0, 150, 25), label='Group 4: med = {}'.format(np.median(star4_med_disp_R)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
# plt.legend(loc=2, frameon=False)
# plt.xlabel(r'$\sigma_{R}$')
# plt.savefig('/Users/amandaquirk/Desktop/med_r_hist.png', bbox_inches='tight')
# plt.close()

#histograms of anisotropies
single_plot()
plt.hist(z_r_1, bins=np.linspace(0, 2, 25), label='Group 1: med = {}'.format(np.median(z_r_1)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(z_r_2, bins=np.linspace(0, 2, 25), label='Group 2: med = {}'.format(np.median(z_r_2)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
plt.hist(z_r_3, bins=np.linspace(0, 2, 25), label='Group 3: med = {}'.format(np.median(z_r_3)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
plt.hist(z_r_4, bins=np.linspace(0, 2, 25), label='Group 4: med = {}'.format(np.median(z_r_4)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
plt.legend(loc=2, frameon=False)
plt.xlabel(r'$\frac{\sigma_{z}}{\sigma_{r}} $')
plt.savefig('/Users/amandaquirk/Desktop/med_z_r_hist.png', bbox_inches='tight')
plt.close()

single_plot()
plt.hist(phi_r_1, bins=np.linspace(0, 2, 25), label='Group 1: med = {}'.format(np.median(phi_r_1)), normed=1, histtype='step', linewidth=1.6,linestyle='--',stacked=True,fill=False, color='b')
plt.hist(phi_r_2, bins=np.linspace(0, 2, 25), label='Group 2: med = {}'.format(np.median(phi_r_2)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=False, color='m', hatch='//', alpha=0.4)
plt.hist(phi_r_3, bins=np.linspace(0, 2, 25), label='Group 3: med = {}'.format(np.median(phi_r_3)), normed=1, histtype='step', linewidth=1.6,stacked=True,fill=True, color='green', alpha=0.4)
plt.hist(phi_r_4, bins=np.linspace(0, 2, 25), label='Group 4: med = {}'.format(np.median(phi_r_4)), normed=1, histtype='step', linewidth=1.6,linestyle='-',stacked=True,fill=False, color='r')
plt.legend(loc=2, frameon=False)
plt.xlabel(r'$\frac{\sigma_{\phi}}{\sigma_{R}}$')
plt.savefig('/Users/amandaquirk/Desktop/med_phi_r_hist.png', bbox_inches='tight')
plt.close()


