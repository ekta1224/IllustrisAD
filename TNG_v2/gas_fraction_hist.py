import numpy as np 
import matplotlib.pyplot as plt

halos = [455957, 458864, 466746, 467575, 478298, 479264, 487042, 489593, 499113, 504142]
for halo in halos:
	print(halo)
	#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
	gas_fraction = np.loadtxt('/Volumes/Titan/analogs/TNGdata/rotated/modified_rotation/{}_gas_properties_TNGv2_rotated.txt'.format(int(halo)), usecols=(7), unpack = True)
	print(np.median(gas_fraction))
	plt.hist(gas_fraction, bins=np.linspace(0.3, 1, 20), alpha=0.5)

plt.xlabel('Neutral Gas Fraction')
plt.savefig('/Users/amandaquirk/Desktop/TNGv2_nH_frac_zoom.png')