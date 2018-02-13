from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import genfromtxt


a=genfromtxt('/Users/dmelgar/Amatrice2016/afters/afters.txt',delimiter='|')
f1=genfromtxt(u'/Users/dmelgar/Slip_inv/Amatrice_NP1/data/model_info/NP1.fault')
f2=genfromtxt(u'/Users/dmelgar/Slip_inv/Amatrice_NP2/data/model_info/NP2.fault')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(f1[:,1], f1[:,2],-f1[:,3],c='r',s=80)
ax.scatter(f2[:,1], f2[:,2],-f2[:,3],c='b',s=80)
ax.legend(['West dip','East dip'])
ax.scatter(a[:,3], a[:,2],-a[:,4],c='k',s=40)

plt.show()