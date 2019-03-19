#still investigating odd LOS velocity patterns for some of the M31 analogs
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

weirdids = [356427, 365141, 383739, 385230, 385343, 387913, 388958, 389486, 390669, 392084, 394566, 395868, 396624, 397062, 397148, 397204, 397889, 398622, 400733, 401663]
goodids = [348189, 361428, 368428, 383612, 388303, 392544, 394155, 394942, 396789, 400179, 401705, 402841, 404964, 405300] # controls

for id in goodids:
    plt.figure(figsize=(12,6))

    gas = np.loadtxt('../data/M31analog_%s_gas_properties.txt'%id)
    x = gas[:,0] 
    y = gas[:,1] 
    z = gas[:,2] 

    stars = np.loadtxt('../data/M31analog_%s_star_properties.txt'%id)
    x2= stars[:,0] 
    y2 = stars[:,1] 
    z2 = stars[:,2] 


    plt.subplot(121)
    plt.hist2d(x, y, bins=100, cmap='viridis',norm=matplotlib.colors.LogNorm())
    cb = plt.colorbar()
    cb.set_label('counts in bin')
    plt.title('gas')

    plt.subplot(122)
    plt.hist2d(x2, y2, bins=100, cmap='magma',norm=matplotlib.colors.LogNorm())
    cb = plt.colorbar()
    cb.set_label('counts in bin')
    plt.title('stars')
    plt.savefig('%s_3d_particles_control.png'%int(id), dpi=300)
    #plt.show()

   #  ax.scatter(x,y,z, c='r', marker='o')
#     ax.scatter(x2,y2,z2, c='b', marker='^')
#     ax.set_xlabel('X Label')
#     ax.set_ylabel('Y Label')
#     ax.set_zlabel('Z Label')

    plt.tight_layout()
    plt.close()
    
