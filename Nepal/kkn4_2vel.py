from obspy import read
from numpy import r_,diff,mean
from datetime import timedelta
from matplotlib import pyplot as plt

fcorner=0.5
Td=timedelta(seconds=45)

sta='SYBC'

n=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXN.sac')
e=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXE.sac')
u=read(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.LXZ.sac')

def low_pass_filter(tr,fcorner,dt,order):
    from scipy.signal import filtfilt,butter
    fnyquist=1./(2*dt)
    print fnyquist
    Fc=fcorner/fnyquist
    print Fc
    b, a = butter(order, Fc,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y
    
n[0].data=low_pass_filter(n[0].data,fcorner,n[0].stats.delta,2)
e[0].data=low_pass_filter(e[0].data,fcorner,e[0].stats.delta,2)
u[0].data=low_pass_filter(u[0].data,fcorner,u[0].stats.delta,2)

n[0].data=n[0].data-mean(n[0].data)
e[0].data=e[0].data-mean(e[0].data)
u[0].data=u[0].data-mean(u[0].data)

n[0].data=r_[0,diff(n[0].data)/n[0].stats.delta]
e[0].data=r_[0,diff(e[0].data)/e[0].stats.delta]
u[0].data=r_[0,diff(u[0].data)/u[0].stats.delta]

n[0].trim(starttime=n[0].stats.starttime,endtime=n[0].stats.starttime+Td)
e[0].trim(starttime=e[0].stats.starttime,endtime=e[0].stats.starttime+Td)
u[0].trim(starttime=u[0].stats.starttime,endtime=u[0].stats.starttime+Td)

n.write(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.vel.n',format='SAC')
e.write(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.vel.e',format='SAC')
u.write(u'/Users/dmelgar/Nepal2015/GPS/cut/'+sta+'.vel.u',format='SAC')

plt.figure()
plt.subplot(311)
plt.plot(n[0].times(),n[0].data)
plt.ylabel('North')
plt.subplot(312)
plt.plot(e[0].times(),e[0].data)
plt.ylabel('East')
plt.subplot(313)
plt.plot(u[0].times(),u[0].data)
plt.ylabel('Up')
plt.show()