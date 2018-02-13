from numpy import genfromtxt
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np

f=genfromtxt(u'/Users/dmelgar/Tohoku2011/Minsons/rupt.txt')

npts = 200
x = f[:,2]
y = f[:,1]
z = (f[:,8]**2+f[:,9]**2)**0.5
# define grid.
xi = np.linspace(x.min(), x.max(), 1000)
yi = np.linspace(y.min(), y.max(), 1000)
# grid the data.
zi = griddata(x, y, z, xi, yi, interp='linear')
# contour the gridded data, plotting dots at the nonuniform data points.
CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
CS = plt.contourf(xi, yi, zi, 15,cmap=plt.cm.rainbow,
                  vmax=70, vmin=0)
                  
plt.colorbar()  # draw colorbar
CS = plt.contour(xi, yi, zi, levels=[10], linewidths=2.0, colors='k')                
                  

# plot data points.
#plt.scatter(x, y, marker='o', s=5, zorder=10)
plt.title('griddata test (%d points)' % npts)
plt.show()