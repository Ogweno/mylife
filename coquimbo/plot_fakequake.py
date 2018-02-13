from numpy import genfromtxt,sqrt
from matplotlib import pyplot as plt

f1=genfromtxt(u'/Volumes/Illapel/FQ/illapel/output/ruptures/illapel.000007.rupt')
f2=genfromtxt(u'/Users/dmelgar/Coquimbo2015/coquimbo_FQ_3km.rupt')

plt.figure()
plt.subplot(121)
plt. scatter(f2[:,1],f2[:,2],c=sqrt(f2[:,8]**2+f2[:,9]**2),lw=0,s=5,cmap=plt.cm.nipy_spectral_r,vmax=16.5)
plt.colorbar()
plt.subplot(122)
plt. scatter(f1[:,1],f1[:,2],c=sqrt(f1[:,8]**2+f1[:,9]**2),lw=0,s=5,cmap=plt.cm.nipy_spectral_r,vmax=16.5)
plt.colorbar()
plt.show()