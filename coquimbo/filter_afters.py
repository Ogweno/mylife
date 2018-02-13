from numpy import genfromtxt,savetxt,c_
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


fault=genfromtxt('/Users/dmelgar/Coquimbo2015/coquimbo_slab1.0.txt')
xs=fault[:,0]
ys=fault[:,1]
zs=fault[:,2]

#Get afters
afters=genfromtxt(u'/Users/dmelgar/Coquimbo2015/afters/afters_3weeks.txt',usecols=[3,4,5])
xa=-afters[:,1]
ya=-afters[:,0]
za=-afters[:,2]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs)
ax.scatter(xa, ya, za,c='r')

ax.set_xlabel('Lon')
ax.set_ylabel('Lat')
ax.set_zlabel('Depth (km)')

plt.show()

savetxt(u'/Users/dmelgar/code/GMT/coquimbo/afters.txt',c_[xa,ya,za],fmt='%.4f\t%.4f\t%.1f')

