from numpy import genfromtxt,where
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

min_mag=1.0

af1=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af2=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records2.txt',usecols=[0,1,2,3,4,5])
i=where(af1[:,5]>min_mag)[0]
j=where(af2[:,5]>min_mag)[0]

print 'Selceted %d events' % (len(i))
xs1=af1[i,3]
ys1=af1[i,2]
zs1=af1[i,4]
xs2=af2[j,3]
ys2=af2[j,2]
zs2=af2[j,4]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs1, ys1, -zs1,s=5,lw=0,c='b')
ax.scatter(xs2, ys2,-zs2,s=5,lw=0,c='b')

ax.set_xlabel('Lon')
ax.set_ylabel('Lat')
ax.set_zlabel('Depth (km)')

plt.show()