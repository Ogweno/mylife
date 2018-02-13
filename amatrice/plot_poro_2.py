from numpy import genfromtxt
from matplotlib import pyplot as plt

c1=genfromtxt('/Users/dmelgar/code/pegrn/output/norcia/all008.dat')
c2=genfromtxt('/Users/dmelgar/code/pegrn/output/norcia_deep/all008.dat')


cmax=0.005
cmin=-0.005
size=30
column=18

plt.figure()

plt.subplot(131)
plt.scatter(c1[:,1],c1[:,0],s=size,c=-c1[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('Normal')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(132)
plt.scatter(c2[:,1],c2[:,0],s=size,c=-c2[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
plt.title('Deep')
cb=plt.colorbar()
cb.set_label('MPa')

plt.subplot(133)
plt.scatter(c2[:,1],c2[:,0],s=size,c=-c2[:,column]/1e6+c1[:,column]/1e6,cmap=plt.cm.seismic,lw=0,vmin=cmin,vmax=cmax)
cb=plt.colorbar()
plt.title('Deep - normal')
cb.set_label('MPa')



plt.show()