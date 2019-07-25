import numpy as np
import os
from os import path

ids = np.loadtxt('M31analogs_halo_props_strictly_noM33.txt')
print(len(ids))
need_data = []
have_data = []
for i in ids:
    res = (str(path.exists('M31analog_%s_gas_properties.txt'%int(i))))
    print (res)
    if res =='False':
        need_data.append(i)
    else:
         have_data.append(i)
#         os.system('git add M31analog_%s_gas_properties.txt'%int(i))
#         os.system('git add M31analog_%s_star_properties.txt'%int(i))
print(need_data)
print(len(need_data))

np.savetxt('get_pdata_for_these_M31analogs.txt', np.column_stack((need_data)), delimiter="  ")

np.savetxt('have_pdata_for_these_M31analogs.txt', np.column_stack((have_data)), delimiter="  ")
#print str(path.exists('M31analog_0_gas_properties.txt'))
