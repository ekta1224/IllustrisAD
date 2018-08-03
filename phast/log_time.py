import numpy as np
from cosmo_tools import snapnum2z, time

ts = np.arange(8.6, 10.4, 0.2)
print ts
ts = ts[::-1]
print ts
for t in ts:
    print t
    print (10.**(t) - 10.**(t-0.2))/1e9


snaps = np.arange(126, 135, 1)
for snap in snaps:
    print snap
    print time(snapnum2z(snap))

#print time(snapnum2z(121))
