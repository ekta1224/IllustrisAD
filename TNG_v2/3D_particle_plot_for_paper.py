from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    levels = [one_sigma, two_sigma][::-1]
    labels = [0.68, 0.95][::-1]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    fmt = {}
    strs = ['0.68', '0.95'][::-1]

    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels[1], origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        #ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    return contour

id = 449972

plt.figure(figsize=(10,4))

gas = np.loadtxt('/Users/ektapatel/Desktop/rotated_data/%s_gas_properties_TNGv2_rotated.txt'%id)
nh = gas[:,7]
mask = np.where(nh > 0.6)
x = gas[:,0][mask]
y = gas[:,1][mask] 
z = gas[:,2][mask] 
print(len(x))

stars = np.loadtxt('/Users/ektapatel/Desktop/rotated_data/%s_star_properties_TNGv2_rotated.txt'%id)
x2= stars[:,0] 
y2 = stars[:,1] 
z2 = stars[:,2] 


ax1 = plt.subplot(121)
print(len(x))
circle1 = plt.Circle((0, 0), 20., color='orange', fill=False, lw=1.5)
circle2 = plt.Circle((0, 0), 20., color='white', fill=False, lw=1.5)
rs = np.array([np.sqrt(x**2.+y**2.) for x,y in zip(x[z<=10],y[z<=10])])
print(len(rs), len(rs[rs<20.]), len(rs[rs<20.])/float(len(rs)))
plt.hist2d(x[z<=10.], y[z<=10.], bins=int(np.sqrt(len(x[z<=10]))), cmap='viridis',norm=matplotlib.colors.LogNorm())
cb = plt.colorbar()
cb.set_label('cells/bin')
print(int(np.sqrt(len(x[z<=10.]))),int(np.sqrt(len(y[z<=10.])))) #150, 150
#density_contour(np.array(x[z<=10.]),np.array(y[z<=10.]),150,150,ax=ax1, colors=['white'], linewidths=[0.75])
ax1.add_artist(circle1)
plt.title('gas')
plt.xlim(-100,100)
plt.ylim(-100,100)
plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')

print(len(np.histogram_bin_edges(x2,bins='stone')))
print(len(x2))
ax2 = plt.subplot(122)
plt.hist2d(x2[z2<=10.], y2[z2<=10.], bins=int(np.sqrt(len(x2[z2<=10.]))), cmap='magma',norm=matplotlib.colors.LogNorm())
#density_contour(x2,y2,int(np.sqrt(len(x2[z2<=10.]))), int(np.sqrt(len(x2[z2<=10.]))),ax=ax2, colors=['white'], linewidths=[0.75])
rs2 = np.array([np.sqrt(x**2.+y**2.) for x,y in zip(x2[z2<=10],y2[z2<=10])])
print(len(rs2), len(rs2[rs2<20.]), len(rs2[rs2<20.])/float(len(rs2)))
ax2.add_artist(circle2)
cb = plt.colorbar()
cb.set_label('particles/bin')
plt.title('stars')

plt.xlim(-100,100)
plt.ylim(-100,100)

plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')


plt.tight_layout()
plt.savefig('%s_3d_rotated_TNGv2.png'%int(id), dpi=300)
plt.close()

