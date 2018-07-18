#columns for data files

#gas file
column 0,1,2: x,y,z position in ckpc/h (divide by h to get to kpc)
column 3,4,5: x,y,z position in km sqrt(a)/s (multiply by sqrt(a) to get to km/s -- for z=0 snapshot this doesn't matter)
column 6: mass in 10^10 Msun/h (multiply by 10^10 and divide by h to get to Msun)
column 7: neutral hydrogen fraction (unitless)

#stars file
column 0,1,2: x,y,z position in ckpc/h (divide by h to get to kpc)
column 3,4,5: x,y,z position in km sqrt(a)/s (multiply by sqrt(a) to get to km/s -- for z=0 snapshot this doesn't matter)
column 6: mass in 10^10 Msun/h (multiply by 10^10 and divide by h to get to Msun)
column 7: stellar age (time of formation in units of the scale factor - use function in redshift2time to convert from z to lookback time)

#dark matter file
column 0,1,2: x,y,z position in ckpc/h (divide by h to get to kpc)
column 3,4,5: x,y,z position in km sqrt(a)/s (multiply by sqrt(a) to get to km/s -- for z=0 snapshot this doesn't matter)
column 6: mass in 10^10 Msun/h (multiply by 10^10 and divide by h to get to Msun)


#Notes:
-all positions should be subtracted from the center of mass (COM) position of the halo itself
-h = 0.704
