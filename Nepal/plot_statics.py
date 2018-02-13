from matplotlib import pyplot as plt
from numpy import genfromtxt,zeros
from glob import glob


sta=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_ALOS/data/station_info/nepal_lrgps_t048.gflist',usecols=0,dtype='S')
lon=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_ALOS/data/station_info/nepal_lrgps_t048.gflist',usecols=[1],dtype='f')
lat=genfromtxt('/Users/dmelgar/Slip_inv/Nepal_ALOS/data/station_info/nepal_lrgps_t048.gflist',usecols=[2],dtype='f')
n=zeros(11)
e=zeros(11)

plt.figure()
for k in range(11):
    print sta[k]
    neu=genfromtxt(u'/Users/dmelgar/Slip_inv/Nepal_ALOS/data/statics/'+sta[k]+'.neu')
    n[k]=neu[0]
    e[k]=neu[1]
    plt.annotate(sta[k],xy=(lon[k],lat[k]))
plt.scatter(lon[0:11],lat[0:11])
plt.quiver(lon[0:11],lat[0:11],e,n,color='green',scale=5,width=0.005)
plt.show()