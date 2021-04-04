import requests 
import numpy as np 
import h5py 

#function that gives path to data
def get(path, params=None):
	# make HTTP GET request to path
	headers = {"api-key":"624abff74a436f05f87b576b68b19a76"}
	r = requests.get(path, params=params, headers=headers)
	
	# raise exception if response code is not HTTP SUCCESS (200)
	r.raise_for_status()
	
	if r.headers['content-type'] == 'application/json':
	    return r.json() # parse json responses automatically
	
	if 'content-disposition' in r.headers:
	    filename = r.headers['content-disposition'].split("filename=")[1]
	    with open(filename, 'wb') as f:
	        f.write(r.content)
	    return filename # return the filename string
	
	return r


#function to get merger trees to deliver snaps and ids as a function of time
def get_ids_snaps(subhalo_id, snap_num):
        url = "http://www.illustris-project.org/api/Illustris-1/snapshots/" + str(snap_num) + "/subhalos/" + str(int(subhalo_id))
        print(url)
        sub = get(url)

        #access merger trees
        mpb1 = get(sub['trees']['sublink_mpb'])
        t = h5py.File(mpb1,'r')
        print(t.keys())
        snapnums = t['SnapNum']
        subfindIDs = t['SubfindID']
        print(snapnums,subfindIDs)

        mysnaps = [134., 133., 132., 131., 129., 126., 121., 112., 98.]
        myids = subfindIDs[np.in1d(snapnums, mysnaps)]
        return myids, mysnaps

#function to get the data I want
def save_particle_data(masterid, subhalo_id, snap_num):
	#para keys found: https://www.illustris-project.org/data/docs/specifications/#sec2a
        params = {'stars':'Coordinates,Velocities,GFM_StellarFormationTime,Masses,GFM_Metallicity','gas':'Coordinates,Velocities,Masses,NeutralHydrogenAbundance,StarFormationRate,GFM_Metallicity','dm':'Coordinates,Velocities,Potential'}

        #url = "http://www.illustris-project.org/api/Illustris-1/snapshots/{}/subhalos/".format(snap_num) + str(int(subhalo_id))
        url = "http://www.illustris-project.org/api/Illustris-1/snapshots/" + str(snap_num) + "/subhalos/" + str(int(subhalo_id))
        print(url)
        sub = get(url)
        #this automatically saves the cutout to your folder, but I don't know how to change that
        saved_filename = get(url + "/cutout.hdf5",params) #hdf5 file with all the data
        f = h5py.File(saved_filename, 'r') #read the file in order to save the data
        
        stars = f['PartsType4']
        gas = f['PartType0']
        dm = f['PartType1']
        subhalo_id = int(subhalo_id)
        
        
        #star_x = f['PartType4']['Coordinates'][:,0]
        # star_y = f['PartType4']['Coordinates'][:,1]
        # star_z = f['PartType4']['Coordinates'][:,2]
        # star_vx = f['PartType4']['Velocities'][:,0]
        # star_vy = f['PartType4']['Velocities'][:,1]
        # star_vz = f['PartType4']['Velocities'][:,2]
        # star_mass = f['PartType4']['Masses']
        # star_tform = f['PartType4']['GFM_StellarFormationTime']
        # star_met = f['PartType4']['GFM_Metallicity']

        gas_x = f['PartType0']['Coordinates'][:,0]
        gas_y = f['PartType0']['Coordinates'][:,1]
        gas_z = f['PartType0']['Coordinates'][:,2]
        gas_vx = f['PartType0']['Velocities'][:,0]
        gas_vy = f['PartType0']['Velocities'][:,1]
        gas_vz = f['PartType0']['Velocities'][:,2]
        gas_mass = f['PartType0']['Masses']
        gas_nf = f['PartType0']['NeutralHydrogenAbundance']
        gas_sfr = f['PartType0']['StarFormationRate']
        gas_met = f['PartType0']['GFM_Metallicity']

        # dm_x = f['PartType1']['Coordinates'][:,0]
        # dm_y = f['PartType1']['Coordinates'][:,1]
        # dm_z = f['PartType1']['Coordinates'][:,2]
        # dm_vx = f['PartType1']['Velocities'][:,0]
        # dm_vy = f['PartType1']['Velocities'][:,1]
        # dm_vz = f['PartType1']['Velocities'][:,2]
        # dm_mass = np.ones(len(dm_vz)) * 6.3 * 10**6

        #save the data to text files -- sorry Bruno 
        #np.savetxt('{}_star.txt'.format(int(subhalo_id)), np.c_[star_x, star_y, star_z, star_vx, star_vy, star_vz, star_mass, star_tform, star_met])
        np.savetxt('./SFHs/M31analog_{}_gas_properties_snap{}.txt'.format(int(masterid), snap_num), np.c_[gas_x, gas_y, gas_z, gas_vx, gas_vy, gas_vz, gas_mass, gas_nf, gas_sfr, gas_met])
        #np.savetxt('{}_dm.txt'.format(int(subhalo_id)), np.c_[dm_x, dm_y, dm_z, dm_vx, dm_vy, dm_vz, dm_mass])
        return

if __name__ == "__main__":
        halo_ids = np.loadtxt('M31analogs_noMM8Gyr_mstar.txt') #61
        print(len(halo_ids))
        
        #halo_ids = np.loadtxt('M31analogs_noMM8Gyr_mstar_noM33.txt')
        #print(len(halo_ids))
        
        #halo_ids = np.loadtxt('M31analogs_MM1_4Gyr_mstar.txt')
        import time


        for halo in halo_ids:
                allids, allsnaps = get_ids_snaps(halo, 135)
                break
        assert False

        for halo in halo_ids:
                start = time.time()
                allids, allsnaps = get_ids_snaps(halo, 135)
                print(allids, allsnaps)
                for id,snap in zip(allids,allsnaps):
                        save_particle_data(halo, id, int(snap))
                
                stop = time.time()
                print(stop-start)
