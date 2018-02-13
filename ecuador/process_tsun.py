from obspy import read
from datetime import timedelta
from obspy.core import UTCDateTime
from obspy import Stream,Trace
from numpy import genfromtxt,where,arange
from scipy.interpolate import interp1d


time_epi=UTCDateTime('2016-04-16T23:58:36.920')
#st=read(u'/Users/dmelgar/Ecuador2016/tsunami/lali.sac')
#st=read('/Users/dmelgar/Ecuador2016/tsunami/32411.sac')
st=read('/Users/dmelgar/Ecuador2016/tsunami/32413.sac')

fcorner=1./7200

def filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt
    
st2=st.copy()
st2[0].data=filter(st[0].data,fcorner,1./st[0].stats.delta,2)

st2.trim(starttime=time_epi-timedelta(hours=8),endtime=time_epi+timedelta(hours=8))


#Read in gauges
g=genfromtxt('/Users/dmelgar/Tsunamis/Ecuador_noshallow/_output/fort.gauge')
tg=g[:,2]
etag=g[:,6]
i=where(etag>0.5)[0]
etag[i]=0
tg_interp=arange(0,tg.max(),15)
I=interp1d(tg,etag)
eta=I(tg_interp)


st_eta=Stream(Trace())
st_eta[0].data=eta
st_eta[0].stats.starttime=time_epi
st_eta[0].stats.delta=15.
st_eta[0].trim(starttime=time_epi-timedelta(hours=8),endtime=time_epi+timedelta(hours=8),pad=True,fill_value=0)