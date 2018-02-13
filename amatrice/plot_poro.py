from numpy import genfromtxt
from matplotlib import pyplot as plt

c0=genfromtxt('/Users/dmelgar/code/pegrn/output/all000.dat')
c1=genfromtxt('/Users/dmelgar/code/pegrn/output/all001.dat')
c2=genfromtxt('/Users/dmelgar/code/pegrn/output/all002.dat')
c3=genfromtxt('/Users/dmelgar/code/pegrn/output/all003.dat')
c4=genfromtxt('/Users/dmelgar/code/pegrn/output/all004.dat')
c5=genfromtxt('/Users/dmelgar/code/pegrn/output/all005.dat')
c6=genfromtxt('/Users/dmelgar/code/pegrn/output/all006.dat')
c7=genfromtxt('/Users/dmelgar/code/pegrn/output/all007.dat')
c8=genfromtxt('/Users/dmelgar/code/pegrn/output/all008.dat')

cmax=0.1
cmin=-0.1
size=30
column=18

plt.figure()

plt.subplot(331)
plt.scatter(c0[:,1],c0[:,0],s=size,c=-c0[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('0 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(332)
plt.scatter(c1[:,1],c1[:,0],s=size,c=-c1[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('1 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(333)
plt.scatter(c2[:,1],c2[:,0],s=size,c=-c2[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
cb=plt.colorbar()
plt.title('2 mo. after EQ')
cb.set_label('MPa')

plt.subplot(334)
plt.scatter(c3[:,1],c3[:,0],s=size,c=-c3[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('3 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(335)
plt.scatter(c4[:,1],c4[:,0],s=size,c=-c4[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('4 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(336)
plt.scatter(c5[:,1],c5[:,0],s=size,c=-c5[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
cb=plt.colorbar()
cb.set_label('MPa')
plt.title('5 mo. after EQ')

plt.subplot(337)
plt.scatter(c6[:,1],c6[:,0],s=size,c=-c6[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('6 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(338)
plt.scatter(c7[:,1],c7[:,0],s=size,c=-c7[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('9 mo. after EQ')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(339)
plt.scatter(c8[:,1],c8[:,0],s=size,c=-c8[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
cb=plt.colorbar()
cb.set_label('MPa')
plt.title('12 mo. after EQ')

plt.figure()
plt.scatter(c8[:,1],c8[:,0],s=size,c=c0[:,column]/1e6-c8[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
cb=plt.colorbar()
cb.set_label('MPa')
plt.title('Difference, 12 mo. - 0 mo.after EQ')

plt.show()