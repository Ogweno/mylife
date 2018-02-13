from numpy import genfromtxt,zeros,c_,mean,arange,where
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace

plt.close("all")

station='crot'
fin='/Users/dmelgar/Lefkada2015/tsunami/raw/'+station+'_ioc_long.txt'
fout='/Users/dmelgar/Lefkada2015/tsunami/sac/'+station+'_long.sac'
plot_out='/Users/dmelgar/Lefkada2015/tsunami/plots/'+station+'.pdf'
spike_threshold=-0.25
dt=60


time=c_[genfromtxt(fin,usecols=0,dtype='S'),genfromtxt(fin,usecols=1,dtype='S')]
t1=UTCDateTime(time[0,0]+time[0,1])
#Make every time relative ins econds to t1
t=zeros(len(time))
for k in range(len(time)):
    t[k]=UTCDateTime(time[k,0]+time[k,1])-t1
#read data
prs=genfromtxt(fin,usecols=2)
#rad=genfromtxt(fin,usecols=3)
#plt.plot(t,prs-mean(prs),t,rad-mean(rad));plt.legend(['prs','rad']) ;plt.show()
plt.plot(t,prs-mean(prs));plt.legend(['prs']) ;plt.show()
#pick one
data=prs.copy()
#data=rad.copy()
#De-spike
i=where(data<spike_threshold)[0]
data[i]=data[i-1]

#Interpolate to regular itnerval
t_in=arange(t[0],t[-1],60)
f=interp1d(t,data)
data_in=f(t_in)
#Processing
data_in=data_in-mean(data_in)
#Put in sac file
st=Stream(Trace())
st[0].stats.station=station
st[0].stats.delta=dt
st[0].stats.starttime=t1
st[0].data=data_in
st.write(fout,format='SAC')
