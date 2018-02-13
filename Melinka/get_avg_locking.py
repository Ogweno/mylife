from numpy import genfromtxt,linspace, meshgrid,c_,where
from mudpy.view import plot_grd
from matplotlib import pyplot as plt
from scipy.interpolate import griddata

fault=genfromtxt('/Users/dmelgar/Slip_inv/Melinka_usgs/output/inverse_models/models/gsi_vr2.6.0011.inv.total')

grdfile='/Users/dmelgar/code/GMT/Melinka/lock.grd'
Xlock,Ylock,lock=plot_grd(grdfile,[0,1],plt.cm.magma,flip_lon=False,return_data=True)

#Interpolate
x=linspace(-75,-73,100)
y=linspace(-44,-42,100)
X,Y=meshgrid(x,y)
z = griddata(fault[:,1:3], (fault[:,8]**2+fault[:,9]**2)**0.5, (X, Y), method='linear',fill_value=0)

#get 1m contour
plt.contour(x, y, z,levels=[1,2,3,4,5,6],lw=0.5)
cs=plt.contour(x, y, z,levels=[1],lw=10)
plt.xlim([-75,-73])
plt.ylim([-44,-42])



path=cs.collections[0].get_paths()
p=path[1]

points=c_[Xlock.ravel(),Ylock.ravel()]
i=where(p.contains_points(points)==True)[0]

m=lock.ravel()[i].mean()

plt.title('Mean locking inside 1m contour is %.2f' % (m))

plt.figure()
plt.scatter(points[i,0],points[i,1],c=lock.ravel()[i],lw=0,vmin=0,vmax=1.0,cmap=plt.cm.magma)
plt.colorbar()
plt.title('Locking inside 1m contour')
plt.xlim([-75,-73])
plt.ylim([-44,-42])

plt.show()