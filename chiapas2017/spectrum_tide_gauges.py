import matplotlib.pyplot as plt
from obspy import read
from numpy import genfromtxt,arange
from scipy.interpolate import interp1d
from obspy import UTCDateTime,Stream,Trace
from mtspec import mtspec
from datetime import timedelta

dt=60.
time_epi=UTCDateTime('2017-09-08T04:49:21')
trim_time=timedelta(seconds=24*3600)

h1=-8.4
h2=-2
shoal=(h1/h2)**0.25
data=read(u'/Users/dmelgar/Chiapas2017/tsunami/long/pchi.sac')
data[0].trim(starttime=time_epi,endtime=time_epi+trim_time)


syn=genfromtxt(u'/Users/dmelgar/Tsunamis/tehuantepec_gauges/_output/gauge00004.txt')
#resample to regular itnerval
tsyn=arange(syn[0,1],syn[-1,1],dt)
f=interp1d(syn[:,1],syn[:,5])
syn_interp=f(tsyn)*shoal
syn=Stream(Trace())
syn[0].data=syn_interp
syn[0].stats.delta=dt
syn[0].stats.starttime=time_epi
syn[0].trim(starttime=data[0].stats.starttime,pad=True,fill_value=0)

#get spectra
delta=dt/3600.
spec_data, freq_data, jackknife, _, _ = mtspec(
    data=data[0].data, delta=delta, time_bandwidth=3.5,
    number_of_tapers=5, nfft=data[0].stats.npts, statistics=True)
    
spec_syn, freq_syn, jackknife, _, _ = mtspec(
    data=syn[0].data, delta=delta, time_bandwidth=3.5,
    number_of_tapers=5, nfft=syn[0].stats.npts, statistics=True)



ax = plt.subplot(111)
pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width / 2.0, pos1.height / 2.0] 
ax.set_position(pos2) # set a new position