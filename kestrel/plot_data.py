from obspy import read
from numpy import genfromtxt,sin,cos,deg2rad,array,c_
from matplotlib import pyplot as plt

n=read(u'/Users/dmelgar/kestrel/BRIC/BRIC.BK/BYN.00.D/BRIC.BK.BYN.00.D.2016.232')
e=read(u'/Users/dmelgar/kestrel/BRIC/BRIC.BK/BYE.00.D/BRIC.BK.BYE.00.D.2016.232')
z=read(u'/Users/dmelgar/kestrel/BRIC/BRIC.BK/BYZ.00.D/BRIC.BK.BYZ.00.D.2016.232')

n[0].data=n[0].data*100e-6
e[0].data=e[0].data*100e-6
z[0].data=z[0].data*100e-6

yl=[-0.08,0.08]

#sopac
g=genfromtxt('/Users/dmelgar/Downloads/pos_brib_57620_00')
x1=g[:,2]-g[0,2]
y1=g[:,3]-g[0,3]
z1=g[:,4]-g[0,4]
x2=g[:,8]-g[0,8]
y2=g[:,9]-g[0,9]
z2=g[:,10]-g[0,10]

#Rotate to local NEU
lat=deg2rad(37.91940521)
lon=deg2rad(-122.15255493)
R=array([[-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat)],[-sin(lon),cos(lon),0],[cos(lon)*cos(lat),cos(lat)*sin(lon),sin(lat)]])
scripps1=R.dot(c_[x1,y1,z1].T).T
scripps2=R.dot(c_[x2,y2,z2].T).T

plt.subplot(311)
plt.plot(n[0].times(),n[0].data,'k')
plt.plot(scripps1[:,0],c='#1E90FF')
plt.plot(scripps2[:,0],c='#DC143C')
plt.xlim([0,len(y1)])
plt.ylabel('North (m)')
plt.legend(['Kestrel RTX','Scripps 1','Scripps 2'])
plt.ylim(yl)

plt.subplot(312)
plt.plot(e[0].times(),e[0].data,'k')
plt.plot(scripps1[:,1],c='#1E90FF')
plt.plot(scripps2[:,1],c='#DC143C')
plt.xlim([0,len(x1)])
plt.ylabel('East (m)')
plt.ylim(yl)

plt.subplot(313)
plt.plot(z[0].times(),z[0].data,'k')
plt.plot(scripps1[:,2],c='#1E90FF')
plt.plot(scripps2[:,2],c='#DC143C')
plt.xlim([0,len(y1)])
plt.ylabel('Up (m)')
plt.xlabel('Seconds')
plt.ylim(yl)

plt.show()

