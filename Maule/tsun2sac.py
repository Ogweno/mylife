from numpy import genfromtxt,arange,where
from obspy import Stream,Trace,UTCDateTime
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from datetime import timedelta

f=u'/Users/dmelgar/Maule2010/wave_gauges/talc.tsun'
fout="/Users/dmelgar/Slip_inv/maule_gps_tg/data/waveforms/talc_trim.tsun"
time_epi=UTCDateTime('2010-02-27T06:34:14')
dt=60
tmax=7200
tcut=3490 # corr

t=genfromtxt(f,usecols=0)
tsun=genfromtxt(f,usecols=1)

#Interpolate to regular intervals
finterp=interp1d(t,tsun,bounds_error=False)
ti=arange(0,tmax,dt)
tsun_interp=finterp(ti)
i=where(ti>=0)[0]
ti=ti[i]
tsun_interp=tsun_interp[i]
tsun_interp=tsun_interp-tsun_interp[0]

#Put in sac file
tsun_out=Stream(Trace())
tsun_out[0].data=tsun_interp
tsun_out[0].stats.delta=dt
tsun_out[0].stats.starttime=time_epi

#Plot to decide on cutting
plt.figure()
plt.plot(tsun_out[0].times(),tsun_out[0].data)
plt.scatter(tsun_out[0].times(),tsun_out[0].data)
plt.show()

#Cut and save
tsun_out[0].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=tcut))
tsun_out.write(fout,format='SAC')
