from matplotlib import pyplot as plt
from numpy import genfromtxt,zeros
from glob import glob
from matplotlib import cm

sta=genfromtxt('/Users/dmelgar/Slip_inv/Nepal/data/station_info/nepal_alos_t157.sta',usecols=0,dtype='S')
lon=genfromtxt('/Users/dmelgar/Slip_inv/Nepal/data/station_info/nepal_alos_t157.sta',usecols=[1],dtype='f')
lat=genfromtxt('/Users/dmelgar/Slip_inv/Nepal/data/station_info/nepal_alos_t157.sta',usecols=[2],dtype='f')
los=zeros(1577)

plt.figure()
for k in range(len(los)):
    neu=genfromtxt(u'/Users/dmelgar/Slip_inv/Nepal/data/statics/'+sta[k]+'.los')
    los[k]=neu[0]
plt.scatter(lon,lat,c=los,cmap=cm.seismic)
plt.show()