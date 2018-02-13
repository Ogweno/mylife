from numpy import genfromtxt,zeros,c_,mean,arange,where
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace

plt.close("all")

station='kaik'
fin='/Users/dmelgar/NewZealand2016/tsunami/'+station+'.txt'
fout='/Users/dmelgar/NewZealand2016/tsunami/'+station+'.sac'
plot_out='/Users/dmelgar/Lefkada2015/tsunami/plots/'+station+'.pdf'
time_epi=UTCDateTime('2016-11-13 11:02:56')
tcut=timedelta(seconds=6*3600)
tprior=timedelta(seconds=2*3600)
spike_threshold=0.3

def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/fnyquist,'highpass')
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
#rad=genfromtxt(fin,usecols=3)
#plt.plot(t,prs-mean(prs),t,rad-mean(rad));plt.legend(['prs','rad']) ;plt.show()
plt.plot(t,prs-mean(prs));plt.legend(['prs']) ;plt.show()
#pick one
data=prs.copy()
#data=rad.copy()
#De-spike
i=where(abs(data)>spike_threshold)[0]
data[i]=data[i-1]

#Interpolate to regular itnerval
t_in=arange(t[0],t[-1],60)
f=interp1d(t,data)
data_in=f(t_in)
#Processing
#De-spike
data_in=data_in-mean(data_in)
fcorner=0.00005
dt=t_in[2]-t_in[1]
#data_fil=highpass(data_in,fcorner,1/dt,2)
#data_fil=data_in
#Put in sac file
st=Stream(Trace())
st[0].stats.station=station
st[0].stats.delta=dt
st[0].stats.starttime=t1
st[0].data=data_in
st[0].trim(starttime=time_epi-tprior,endtime=time_epi+tcut)
st[0].data=st[0].data-mean(st[0].data[0:20])
st.write(fout,format='SAC')

#Make plot and save
#Plot processed
#plt.close("all")
plt.figure(figsize=(10,3))
plt.plot((st[0].times()-tprior.seconds)/60,st[0].data)
plt.grid()
plt.ylim([-5,5])
plt.xlim([-60,180])
plt.xlabel('Minutes after OT')
plt.ylabel('Sea surface height (m)')
plt.subplots_adjust(bottom=0.2)
plt.title('Station '+station.upper())
plt.savefig(plot_out)
plt.show()