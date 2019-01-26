import numpy as np 
import matplotlib.pyplot as plt 
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 

halo = '361428'
star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass_all, star_factor = np.loadtxt('../data/M31analog_{}_star_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True) 
gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass_all, gas_fraction = np.loadtxt('../data/M31analog_{}_gas_properties_rotated.txt'.format(halo), usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True)

#Illustric value
h = 0.704

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y):#, z):
	#convert to kpc (coordinates are already shifted to 0,0,0)
	#z_kpc = z / h
	return np.sqrt((x**2) + (y**2)) #+ (z_kpc**2)) #km

#calculate the radial velocity
def radial_v(x, y, vx, vy):
	#convert kpc/h to km
	x_km = x * 3.086e+16 / h
	y_km = y * 3.086e+16 / h
	dot_product = (x_km * vx) + (y_km * vy)# + (z * vz) 
	r_mag = distance(x_km, y_km)#, z)
	return dot_product / r_mag #km/s

#calculate the rotation speed, assuming it's the same as tangential once averaged over the radial bin
def vrot(x, y, vx, vy):
	v_rad = radial_v(x, y, vx, vy) #km/s
	v_tot = np.sqrt((vx**2) + (vy**2))# + (vz**2))
	v_rot = np.sqrt((v_tot**2) + (v_rad**2))
	return v_rot #km/s

neutral_gas = gas_fraction > 0.95 #what is the best limit?
gas_x = gas_x[neutral_gas]
gas_y = gas_y[neutral_gas]
gas_vx = gas_vx[neutral_gas]
gas_vy = gas_vy[neutral_gas]

#get rid of star partiles that are really wind
not_wind = star_factor > 0
star_x = star_x[not_wind]
star_y = star_y[not_wind]
star_vx = star_vx[not_wind]
star_vy = star_vy[not_wind]
star_factor = star_factor[not_wind]

#calculate the ages of stars from the scale factor
def star_age(scale_factor):
	zform = [(1./a) -1. for a in scale_factor]
	cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
	tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
	tage = np.array([13.8-t for t in tform])
	return tage

age = star_age(star_factor)

#divide the stars into 4 age bins -- do we want larger bins? -- yes do this later amanda
group1 = (0 <= age) & (age <=.06) #average age 30 Myr 
group2 = (.1 <= age) & (age <=.9) #average age 400 Myr
group3 = (1 <= age) & (age <= 3) #average age 2 Gyr
group4 = (3 <= age) & (age <= 11) #average age 4 Gyr

v = vrot(star_x[group4], star_y[group4], star_vx[group4], star_vy[group4])


plt.scatter(star_x[group4] / h, star_y[group4] / h, c=v, alpha = .5, s=6)#, vmin=100, vmax=350)
clb=plt.colorbar()
clb.set_label('km/s')
plt.xlim(-20,20)
plt.ylim(-20, 20)
plt.xlabel('kpc')
plt.ylabel('kpc')
plt.savefig('/Users/amandaquirk/Desktop/vrot_star_notwind_group4.png')




