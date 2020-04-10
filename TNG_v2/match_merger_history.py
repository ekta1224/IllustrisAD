import h5py 
import numpy as np 

#read in analog IDs -- this matches the index in the merger file
revised_ids = np.loadtxt('/Volumes/Titan/analogs/IllustrisAD/TNG_v2/M31analogs_halo_props_TNG100_revised.txt', usecols=(0,), unpack=True)
revised_ids = [int(a) for a in revised_ids]

#read in the merger file
merger_data = h5py.File('/Volumes/Titan/analogs/TNGdata/merger_history_099.hdf5', 'r')

#pull out the params we want: time of last major merger, number of major mergers, and number of minor mergers and match the data to our analogs
last_major_time = merger_data['SnapNumLastMajorMerger'][revised_ids]
num_major = merger_data['NumMajorMergersTotal'][revised_ids]
num_minor = merger_data['NumMinorMergersTotal'][revised_ids]

#save file
np.savetxt('/Volumes/Titan/analogs/TNGdata/M31analogs_merger_props_TNG100_revised.txt', np.c_[revised_ids, last_major_time, num_major, num_minor], header='analog ID, snapshot number of last major merger, number of total major mergers, number of total minor mergers')

