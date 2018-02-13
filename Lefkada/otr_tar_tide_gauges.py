from numpy import genfromtxt,zeros,c_,mean,arange,where,isnan
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace

plt.close("all")

station='tara'
col=5#2
fin='/Users/dmelgar/Lefkada2015/tsunami/radar 16nov_18nov_2015.TXT'
fout='/Users/dmelgar/Lefkada2015/tsunami/sac/'+station+'.sac'
plot_out='/Users/dmelgar/Lefkada2015/tsunami/plots/'+station+'.pdf'
time_epi=UTCDateTime('2015-11-17T07:10:07')
tcut=timedelta(seconds=6*3600)
tprior=timedelta(seconds=1*3600)
spike_threshold=0.3

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



time=c_[genfromtxt(fin,usecols=0,delimiter=';',dtype='S'),genfromtxt(fin,usecols=1,dtype='S',delimiter=';')]
t1=UTCDateTime(time[0,0].split('-')[2].strip()+'-'+time[0,0].split('-')[1]+'-'+time[0,0].split('-')[0]+'T'+time[0,1][0:2]+':'+time[0,1][3:5])
#Make every time relative ins econds to t1
t=zeros(len(time))
for k in range(len(time)):
    ttemp=UTCDateTime(time[k,0].split('-')[2].strip()+'-'+time[k,0].split('-')[1]+'-'+time[k,0].split('-')[0]+'T'+time[k,1][0:2]+':'+time[k,1][3:5])
    print ttemp
    t[k]=ttemp-time_epi
#read data
prs=genfromtxt(fin,usecols=col,delimiter=';')/100
#rad=genfromtxt(fin,usecols=3)
#plt.plot(t,prs-mean(prs),t,rad-mean(rad));plt.legend(['prs','rad']) ;plt.show()
plt.plot(t,prs-mean(prs));plt.legend(['prs']) ;plt.show()
#pick one
data=prs.copy()
i=where(isnan(data)==1)[0]
data[i]=data[i-1]
#data=rad.copy()
#De-spike
#i=where(abs(data)>spike_threshold)[0]
#data[i]=data[i-1]

#Interpolate to regular itnerval
t_in=arange(t[0],t[-1],60)
f=interp1d(t,data)
data_in=f(t_in)
#Processing
#De-spike
data_in=data_in-mean(data_in)
fcorner=0.00005
dt=t_in[2]-t_in[1]
data_fil=highpass(data_in,fcorner,1/dt,2)
#data_fil=data_in
#Put in sac file
st=Stream(Trace())
st[0].stats.station=station
st[0].stats.delta=dt
st[0].stats.starttime=t1
st[0].data=data_fil
st[0].trim(starttime=time_epi-tprior,endtime=time_epi+tcut)
st[0].data=st[0].data-mean(st[0].data[0:20])
st.write(fout,format='SAC')

#Make plot and save
#Plot processed
#plt.close("all")
plt.figure(figsize=(10,3))
plt.plot((st[0].times()-tprior.seconds)/60,st[0].data)
plt.grid()
plt.ylim([-0.12,0.12])
plt.xlim([-60,180])
plt.xlabel('Minutes after OT')
plt.ylabel('Sea surface height (m)')
plt.subplots_adjust(bottom=0.2)
plt.title('Station '+station.upper())
plt.savefig(plot_out)
plt.show()