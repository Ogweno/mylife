from numpy import genfromtxt,zeros,c_,mean,arange
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace


fin='/Users/dmelgar/Maule2010/wave_gauges/ancu_ioc.txt'
fout='/Users/dmelgar/Maule2010/wave_gauges/ancu.tsun'
time_epi=UTCDateTime('2010-02-27T6:34:14')
tcut=timedelta(seconds=7200)
station='ancu'


def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt



time=c_[genfromtxt(fin,usecols=0,dtype='S'),genfromtxt(fin,usecols=1,dtype='S')]
t1=UTCDateTime(time[0,0]+time[0,1])
#Make every time relative ins econds to t1
t=zeros(len(time))
for k in range(len(time)):
    t[k]=UTCDateTime(time[k,0]+time[k,1])-time_epi
#read data
prs=genfromtxt(fin,usecols=2)
rad=genfromtxt(fin,usecols=2)
plt.plot(t,prs-mean(prs),t,rad-mean(rad));plt.legend(['prs','rad']) ;plt.show()
#pick one
data=prs.copy()
#Interpolate to regular itnerval
t_in=arange(t[0],t[-1],60)
f=interp1d(t,data)
data_in=f(t_in)
#Processing
data_in=data_in-mean(data_in)
fcorner=0.0001
dt=t_in[2]-t_in[1]
data_fil=highpass(data_in,fcorner,1/dt,2)
#Put in sac file
st=Stream(Trace())
st[0].stats.station=station
st[0].stats.delta=dt
st[0].stats.starttime=t1
st[0].data=data_fil
st[0].trim(starttime=time_epi,endtime=time_epi+tcut)
st[0].data=st[0].data-st[0].data[0]
st.write(fout,format='SAC')