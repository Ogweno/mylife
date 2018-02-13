from obspy import read
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.signal import butter,filtfilt
from datetime import timedelta
from numpy import zeros,r_

td=5
fcorner=1./20
e=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.e')
e[0].data=r_[zeros(td),e[0].data]
n=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.n')
n[0].data=r_[zeros(td),n[0].data]
u=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.u')
u[0].data=r_[zeros(td),u[0].data]
ste=read(u'/Users/dmelgar/Slip_inv/Nepal/output/forward_models/xcorr.KATNP.vel.e')
stn=read(u'/Users/dmelgar/Slip_inv/Nepal/output/forward_models/xcorr.KATNP.vel.n')
stu=read(u'/Users/dmelgar/Slip_inv/Nepal/output/forward_models/xcorr.KATNP.vel.u')
time_epi=UTCDateTime('2015-04-25T06:11:26')

def hfilter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y
fnyquist=1./(2*e[0].stats.delta)
Fc=fcorner/fnyquist
ste[0].data=hfilter(ste[0].data,Fc,3)
stn[0].data=hfilter(stn[0].data,Fc,3)
stu[0].data=hfilter(stu[0].data,Fc,3)

e.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))
n.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))
u.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))
ste.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))
stn.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))
stu.trim(starttime=time_epi+timedelta(seconds=5),endtime=time_epi+timedelta(seconds=100))

plt.figure()
plt.subplot(311)
plt.plot(e[0].times(),e[0].data)
plt.plot(ste[0].times(),ste[0].data)
plt.legend(['Data','Synthetic'])
plt.subplot(312)
plt.plot(n[0].times(),n[0].data)
plt.plot(stn[0].times(),stn[0].data)
plt.subplot(313)
plt.plot(u[0].times(),u[0].data)
plt.plot(stu[0].times(),stu[0].data)
plt.show()

e.write(u'/Users/dmelgar/Slip_inv/Nepal/data/waveforms/KATNP.vel.e',format='SAC')
n.write(u'/Users/dmelgar/Slip_inv/Nepal/data/waveforms/KATNP.vel.n',format='SAC')
u.write(u'/Users/dmelgar/Slip_inv/Nepal/data/waveforms/KATNP.vel.u',format='SAC')


