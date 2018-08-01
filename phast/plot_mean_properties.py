import numpy as np
import matplotlib
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

noMM = np.loadtxt('M31analogs_norecentMM_mean_props.txt')
leftz = noMM[:,0]
rightz = noMM[:,1]
leftsfr = noMM[:,2]
rightsfr = noMM[:,3]

zratio_noMM = np.array(leftz)/np.array(rightz)
sfrratio_noMM = np.array(leftsfr)/np.array(rightsfr)

MM = np.loadtxt('M31analogs_recentMM_mean_props.txt')
leftz = MM[:,0]
rightz = MM[:,1]
leftsfr = MM[:,2]
rightsfr = MM[:,3]
rightsfr = [1e-5 if x==0. else x for x in rightsfr]
zratio_MM = np.array(leftz)/np.array(rightz)
sfrratio_MM = np.array(leftsfr)/np.array(rightsfr)

plt.figure()
min = np.min([np.min(zratio_noMM), np.min(zratio_MM)])
max = np.max([np.max(zratio_noMM), np.max(zratio_MM)])
plt.hist(zratio_noMM, histtype='step', bins=6, range=(min,max), label='no recent MM', lw=2)
plt.hist(zratio_MM, histtype='step', bins=6, range=(min,max),label='recent MM', lw=2)
plt.xlabel('metallicity ratio (L/R)')
plt.ylabel('number of analogs')
plt.legend()
plt.savefig('metallicity_ratio.pdf')


plt.figure()
min = np.min([np.min(sfrratio_noMM), np.min(sfrratio_MM)])
max = np.max([np.max(sfrratio_noMM), np.max(sfrratio_MM)])
plt.hist(sfrratio_noMM, histtype='step', bins = np.logspace(np.log10(min), np.log10(max), 6), label='no recent MM', lw=2)
plt.hist(sfrratio_MM, histtype='step', bins = np.logspace(np.log10(min), np.log10(max), 6), label='recent MM', lw=2)
plt.xscale("log")
plt.legend()
plt.xlabel('SFR ratio (L/R)')
plt.ylabel('number of analogs')
plt.savefig('SFR_ratio.pdf')

