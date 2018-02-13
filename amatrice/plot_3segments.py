from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import genfromtxt,r_,array,savetxt

n=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/north_segment.fault',usecols=[1,2,3])
s=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/south_segment_prune.fault',usecols=[1,2,3])
c=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/connector.fault',usecols=[1,2,3])
fullc=genfromtxt('/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/connector.fault')
out='/Users/dmelgar/Slip_inv/Amatrice_3plane/data/model_info/connector_prune.fault'

k=10
i0=array([3,4,5,6])
i1=array([3,4,5,6])+k
i2=array([3,4,5,6])+2*k
i3=array([2,3,4,5,6])+3*k
i4=array([2,3,4,5,6])+4*k
i5=array([2,3,4,5])+5*k
i6=array([2,3,4,5])+6*k
i7=array([1,2,3,4,5])+7*k
i8=array([1,2,3,4,5])+8*k
i9=array([1,2,3,4])+9*k
i10=array([1,2,3,4])+10*k
i11=array([0,1,2,3,4])+11*k

i=r_[i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]
c=c[i,:]
fullc=fullc[i,:]
savetxt(out,fullc,fmt='%d\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(n[:,0], n[:,1], -n[:,2], c='k')
ax.scatter(s[:,0], s[:,1], -s[:,2], c='k')
ax.scatter(c[:,0], c[:,1], -c[:,2], c='r')

plt.show()



