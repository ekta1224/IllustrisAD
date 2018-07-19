import numpy as np 
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u 

'''
1. fix getting position of main halo
2. do first without radial binning
3. how to divide the gas inot age bins
4. do with radial binning
	should we smooth the velocities?
5. make code modular and merge with Ekta's
'''

#load in the positions, velocities, and masses of the stars and gas; the HI fraction, and the age of stars
star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass_all, star_factor = np.loadtxt('../data/M31analog_361428_star_properties_rotated.txt', usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True) 
gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass_all, gas_fraction = np.loadtxt('../data/M31analog_361428_gas_properties_rotated.txt', usecols=(0, 1, 2, 3, 4, 5, 6, 7,), unpack = True)

vx_sys = 2.780412292480468750e+02 #automate this
vy_sys = -1.491532287597656250e+02
vz_sys = -2.593022155761718750e+02

h = 0.704

#calculate distance from center-- we don't need to deproject the position (yay simulations)
def distance(x, y, z):
	#convert to kpc and shift coordinates to be centered on 0,0,0
	x_kpc = x / h
	y_kpc = y / h
	z_kpc = z / h
	return np.sqrt((x_kpc**2) + (y_kpc**2) + (z_kpc**2)) #kpc

def radial_v(vx, vy, vz):
	return np.sqrt((vx**2) + (vy**2) + (vz**2))

vsys = radial_v(vx_sys, vy_sys, vz_sys)

#get rid of star partiles that are really wind
not_wind = star_factor > 0
star_x = star_x[not_wind]
star_y = star_y[not_wind]
star_z = star_z[not_wind]
star_vx = star_vx[not_wind]
star_vy = star_vy[not_wind]
star_vz = star_vz[not_wind]
star_factor = star_factor[not_wind]

star_v = radial_v(star_vx, star_vy, star_vz)

star_r = distance(star_x, star_y, star_z)
gas_r = distance(gas_x, gas_y, gas_z)

#only want particles/cells that correspond to a high HI gas fraction
neutral_gas = gas_fraction > 0.7
gas_r = gas_r[neutral_gas]
gas_x = gas_x[neutral_gas]
gas_y = gas_y[neutral_gas]
gas_vx = gas_vx[neutral_gas]
gas_vy = gas_vy[neutral_gas]
gas_vz = gas_vz[neutral_gas]
gas_v = radial_v(gas_vx, gas_vy, gas_vz)

def PA(x, y):
	return np.arctan(x / y) #radians -- yes x / y because of how the rotation works

star_PA = PA(star_x, star_y)
gas_PA = PA(gas_x, gas_y)

#calculate the rotation speed, assuming spherical
def vrot(v, vsys, PA):
	v_rot = v - vsys * np.sqrt(1 + np.tan(PA)**2)
	return abs(v_rot)

star_vrot = vrot(star_v, vsys, star_PA)
gas_vrot = vrot(gas_v, vsys, gas_PA)

#calculate the ages of stars
def star_age(scale_factor):
	zform = [(1./a) -1. for a in scale_factor]
	cosmo = FlatLambdaCDM(H0=70.4, Om0=0.2726, Ob0=0.0456)
	tform = [float(cosmo.age(z)/u.Gyr) for z in zform] #time the stars formed in lookback time, not age
	tage = np.array([13.8-t for t in tform])
	return tage

age = star_age(star_factor)

#divide the stars into 4 age bins
group1 = (.0 <= age) & (age <=.06) #average age 30 Myr
group2 = (.34 <= age) & (age <=.46) #average age 400 Myr
group3 = (1.7 <= age) & (age <=2.3) #average age 2 Gyr
group4 = (3.7 <= age) & (age <=4.3) #average age 4 Gyr

star1_r = star_r[group1]
star1_vrot = star_vrot[group1]
star1_PA = star_PA[group1]

star2_r = star_r[group2]
star2_vrot = star_vrot[group2]
star2_PA = star_PA[group2]

star3_r = star_r[group3]
star3_vrot = star_vrot[group3]
star3_PA = star_PA[group3]

star4_r = star_r[group4]
star4_vrot = star_vrot[group4]
star4_PA = star_PA[group4]

#divide stars and gas into radial and PA bins and average the velocities within each bin 
#first divide into radial bins


#calculate the asymmetric drift
def va(v_gas, v_star):
	return v_gas - v_star

#plots-- 4 rotation curves and 1 histogram of AD
plt.scatter(gas_r, gas_vrot, c = 'darkgrey', alpha = 0.4)
plt.scatter(star1_r, star1_vrot, c = 'b', alpha = 0.4)
plt.scatter(star2_r, star2_vrot, c = 'm', alpha = 0.4)
plt.scatter(star3_r, star3_vrot, c = 'k', alpha = 0.4)
plt.scatter(star4_r, star4_vrot, c = 'r', alpha = 0.4)
plt.ylim(0,300)
plt.xlim(0,20)
plt.ylabel('Rotation Velocity')
plt.xlabel('r')
plt.show()



