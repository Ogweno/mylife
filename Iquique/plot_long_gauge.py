from numpy import genfromtxt,c_,zeros,mean
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams.update({'font.size': 28})

tsun=genfromtxt('/Users/dmelgar/Iquique2014/tsunami/iquique_long.txt',usecols=2)
tsun=tsun-mean(tsun)
time=c_[genfromtxt('/Users/dmelgar/Iquique2014/tsunami/iquique_long.txt',usecols=0,dtype='S'),genfromtxt('/Users/dmelgar/Iquique2014/tsunami/iquique_long.txt',usecols=1,dtype='S')]
time_epi=UTCDateTime('2014-04-01T23:46:47')
#Make every time relative ins econds to t1
t=zeros(len(time))
for k in range(len(time)):
    t[k]=UTCDateTime(time[k,0]+time[k,1])-time_epi
t=t/3600

plt.figure()
plt.plot(t,tsun)
plt.grid()
plt.xlim([-5,15])
plt.ylim([-3,2.5])
plt.xlabel('Hours after OT',fontsize=28)
plt.ylabel('Sea Surface Height (m)',fontsize=28)
plt.show()
