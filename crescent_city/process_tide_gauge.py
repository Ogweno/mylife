from numpy import genfromtxt,zeros,c_,mean,r_,diff
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from datetime import timedelta
from obspy import Stream,Trace

plt.close("all")

station='cres'
#fin='/Users/dmelgar/tidegauge_noise/crescent_city/cres_01_2017.csv'
#fout='/Users/dmelgar/tidegauge_noise/crescent_city/cres_01_2017.sac'

fin='/Users/dmelgar/tidegauge_noise/crescent_city/cres_09_2015.csv'
fout='/Users/dmelgar/tidegauge_noise/crescent_city/cres_09_2015.sac'

#fin='/Users/dmelgar/tidegauge_noise/eureka/eure_07_2017.csv'
#fout='/Users/dmelgar/tidegauge_noise/eureka/eure_07_2017.sac'

#time_epi=UTCDateTime('2017-01-01T00:00:00')
time_epi=UTCDateTime('2015-09-01T00:00:00')

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



time=c_[genfromtxt(fin,usecols=0,dtype='S',delimiter=',')]
t1=UTCDateTime(time[0,0])
#Make every time relative ins econds to t1
t=zeros(len(time))
for k in range(len(time)):
    t[k]=UTCDateTime(time[k,0])-time_epi
#read data
prs=genfromtxt(fin,usecols=1,delimiter=',')

plt.figure()
plt.plot(t,prs)
plt.figure()
plt.plot(t,r_[360,diff(t)])
plt.show()

data=prs.copy()


dt=t[1]-t[0]

#Put in sac file
st=Stream(Trace())
st[0].stats.station=station
st[0].stats.delta=dt
st[0].stats.starttime=t1
st[0].data=data
st[0].data=st[0].data
st.write(fout,format='SAC')
