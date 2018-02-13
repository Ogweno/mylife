from numpy import genfromtxt,zeros,arange,sin,mean,savetxt,c_
from matplotlib import pyplot as plt
from obspy import read
from obspy.core import UTCDateTime

st=read(u'/Users/dmelgar/NewZealand2016/tsunami/kaik.sac')
fout=u'/Users/dmelgar/code/GMT/NewZealand/kaik.txt'

time_epi=UTCDateTime('2016-11-13 11:02:56')

amplitude=(0.93--0.73)/2
period=7300
phase=-0.05

t=st[0].times()
y=amplitude*sin((1./period)*t+phase)

eta=st[0].data

plt.close('all')
plt.figure()
plt.plot(t,eta)
plt.plot(t,y)
plt.plot(t,eta-y)




#cut
st[0].data=eta-y
st[0].trim(starttime=time_epi-600,endtime=time_epi+14400)
plt.figure()
plt.plot(st[0].times(),st[0].data)

dout=st[0].data-st[0].data[0]
tout=st[0].times()-600
tout=tout/60.

savetxt(fout,c_[tout,dout],fmt='%.3f')

plt.show()