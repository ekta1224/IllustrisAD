from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.colors import LogNorm

h = fits.open('broadband_419510.fits')
data = h[14].data

data_tot = []
for i in range(0, 36):
	data_ind = np.array(data[i])
	#plt.imshow(data_ind, cmap='gray')
	#plt.colorbar()
	#plt.savefig('mock_image_{}.png'.format(i))
	#plt.close()
	data_tot =+ data_ind

plt.imshow(data[21], cmap='gray')#, vmin=0, vmax=17)
plt.colorbar()
plt.ylim(50, 200)
plt.xlim(50, 200)
plt.axis('off')
plt.savefig('/Users/amandaquirk/Desktop/mock_image_21_grey.png', bbox='tight')
plt.close()	