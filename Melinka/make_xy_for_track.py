from numpy import savetxt,linspace,c_,genfromtxt
from matplotlib import pyplot as plt

#fout=u'/Users/dmelgar/Melinka2016/etc/slice360.xy'
#x=linspace(-74.5,-73.5,50)
#y=linspace(-43.2,-43.6,50)
#
#savetxt(fout,c_[x+360,y],fmt='%.4f')
#


#Make plots of slice
hayes=genfromtxt(u'/Users/dmelgar/Melinka2016/etc/hayes_slice.xyz')
moreno=genfromtxt(u'/Users/dmelgar/Melinka2016/etc/moreno_slice.xyz')
planar=genfromtxt(u'/Users/dmelgar/Melinka2016/etc/planar_slice.xyz')


plt.figure()
plt.plot(hayes[:,0]-360,hayes[:,2],c='r',lw=2)
plt.plot(moreno[:,0],moreno[:,2],c='g',lw=2)
plt.plot(planar[:,0],-planar[:,2]+5,c='k',lw=2)
plt.legend(['Hayes','Moreno','Planar'])

plt.show()