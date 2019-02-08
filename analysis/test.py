#testing data
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../data/M31analog_348189_star_properties_rotated.txt')
mass = data[:,6]/0.704

substars = 7.521578311920166016/0.704 #subhalo stellar mass, 1e10 Msun

print np.sum(mass) #1e10 Msun units
print substars, substars*1e10/len(mass)


plt.figure()
plt.hist(mass*1e10)
plt.savefig('348189_stars_masses.pdf')

print len(mass)
print (len(mass)*7e5)/1e10
