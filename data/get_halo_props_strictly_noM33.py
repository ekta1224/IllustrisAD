import numpy as np

master = np.loadtxt('M31analogs_halo_props_noM33.txt')

snom33_ids = np.loadtxt('M31analogs_halo_props_strictly_noM33.txt')

inds = np.in1d(master[:,0], snom33_ids)
print(inds)
print(len(master[inds]))

print(master[inds][0])

stellar_data = np.loadtxt('M31analogs_halo_props_strictly_noM33_stellar_assembly.txt')
#print(stellar_data[:,0])
print(len(stellar_data))

print(len(master[inds][:,0]))

np.savetxt('M31analogs_halo_props_strictly_noM33.txt', np.column_stack((master[inds][:,0], master[inds][:,1], master[inds][:,2], master[inds][:,3], master[inds][:,4], master[inds][:,5], master[inds][:,6], master[inds][:,7], master[inds][:,8], master[inds][:,9], master[inds][:,10], stellar_data[:,1], stellar_data[:,2], stellar_data[:,3], stellar_data[:,4], stellar_data[:,5])), delimiter = " ")



print(stellar_data[:,0][-1])
print(master[inds][:,0][-1])
