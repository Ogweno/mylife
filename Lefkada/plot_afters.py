from matplotlib import pyplot as plt
from numpy import genfromtxt,zeros,sin,cos,deg2rad
from mpl_toolkits.mplot3d import Axes3D
from obspy.core.util.geodetics import gps2DistAzimuth

a=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA.txt')
f=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada70/data/model_info/lefkada65.fault')

lon=a[:,5]
lat=a[:,4]
depth=-a[:,6]
hypo=[20.6002,38.6655]
x=zeros(len(lon))
y=zeros(len(lon))
for k in range(len(lon)):
    dist,az,baz=gps2DistAzimuth(lat[k], lon[k], hypo[1], hypo[0])
    x[k]=dist*sin(deg2rad(baz))/1000
    y[k]=dist*cos(deg2rad(baz))/1000


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(lon,lat,depth)
ax.set_xlabel('lon')
ax.set_ylabel('lat')
ax.set_zlabel('depth (km)')
ax.view_init(azim=-90,elev=90.)
ax.scatter(f[:,1],f[:,2],-f[:,3],c='r',s=100)
ax.set_zlim([0,-20])
ax.set_xlim([20,21])
ax.set_ylim([38,39])

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x,y,depth)
#ax.set_xlabel('Easting (km)')
#ax.set_ylabel('Northing (km)')
#ax.set_zlabel('Depth (km)')
#ax.set_zlim([-20,0])
#ax.set_xlim([-50,50])
#ax.set_ylim([-50,50])
#ax.view_init(azim=-90,elev=90.)
#ax.scatter(0,0,-10.7,s=120,c='r')
plt.show()



