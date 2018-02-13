from datetime import timedelta
from matplotlib import pyplot as plt
from obspy import read,UTCDateTime
from butter import highpass
from numpy import mean

def read_ioc(ioc_file,dt,flip=False):
    from obspy import Stream,Trace,UTCDateTime
    from numpy import genfromtxt,zeros,arange
    from scipy.interpolate import interp1d
    
    month_day=genfromtxt(ioc_file,usecols=0,dtype='S')
    hour=genfromtxt(ioc_file,usecols=1,dtype='S')
    month0=float(month_day[0].split('-')[1])
    day0=float(month_day[0].split('-')[2])
    hour0=float(hour[0].split(':')[0])
    minute0=float(hour[0].split(':')[1])
    second0=float(hour[0].split(':')[2])
    t=zeros(len(hour))
    for k in range(len(t)):
        day_c=float(month_day[k].split('-')[2])
        hour_c=float(hour[k].split(':')[0])
        minute_c=float(hour[k].split(':')[1])
        second_c=float(hour[k].split(':')[2])
        t[k]=(day_c-day0)*(24*60*60)+(hour_c-hour0)*3600+(minute_c-minute0)*60+(second_c-second0)
    data=genfromtxt(ioc_file,usecols=2)
    #Inteprolate
    finterp=interp1d(t,data,bounds_error=False)
    ti=arange(t[0],t[-1],dt)
    tsun_interp=finterp(ti)
    st=Stream(Trace())
    st[0].data=tsun_interp
    st[0].stats.delta=dt
    t1=UTCDateTime(month_day[0]+'T'+hour[0])
    st[0].stats.starttime=t1
    return st

time_epi=UTCDateTime('2010-02-27T06:34:14')
tcut=3490 # corr
corr=read_ioc('/Users/dmelgar/Maule2010/wave_gauges/corr_ioc.txt',60)
talc=read_ioc('/Users/dmelgar/Maule2010/wave_gauges/talc_ioc.txt',60)

def read_ioc(ioc_file,dt,flip=False):
    from obspy import Stream,Trace,UTCDateTime
    from numpy import genfromtxt,zeros,arange
    from scipy.interpolate import interp1d
    
    month_day=genfromtxt(ioc_file,usecols=0,dtype='S')
    hour=genfromtxt(ioc_file,usecols=1,dtype='S')
    month0=float(month_day[0].split('-')[0])
    day0=float(month_day[0].split('-')[1])
    hour0=float(hour[0].split(':')[0])
    minute0=float(hour[0].split(':')[1])
    second0=float(hour[0].split(':')[2])
    t=zeros(len(hour))
    for k in range(len(t)):
        day_c=float(month_day[k].split('-')[0])
        hour_c=float(hour[k].split(':')[0])
        minute_c=float(hour[k].split(':')[1])
        second_c=float(hour[k].split(':')[2])
        t[k]=(day_c-day0)*(24*60*60)+(hour_c-hour0)*3600+(minute_c-minute0)*60+(second_c-second0)
    data=genfromtxt(ioc_file,usecols=2)
    #Inteprolate
    finterp=interp1d(t,data,bounds_error=False)
    ti=arange(t[0],t[-1],dt)
    tsun_interp=finterp(ti)
    st=Stream(Trace())
    st[0].data=tsun_interp
    st[0].stats.delta=dt
    t1=UTCDateTime('2010-02-27T00:00:00')
    st[0].stats.starttime=t1
    return st

valp=read_ioc('/Users/dmelgar/Maule2010/wave_gauges/valp_ioc.txt',60)

fcorner=1./7200

valp[0].data=valp[0].data-mean(valp[0].data)
valp[0].data=highpass(valp[0].data,fcorner,valp[0].stats.delta,2)
corr[0].data=corr[0].data-mean(corr[0].data)
corr[0].data=highpass(corr[0].data,fcorner,corr[0].stats.delta,2)
talc[0].data=talc[0].data-mean(talc[0].data)
talc[0].data=highpass(talc[0].data,fcorner,talc[0].stats.delta,2)

tcut=80*60
valp[0].trim(starttime=time_epi-timedelta(seconds=360),endtime=time_epi+timedelta(seconds=tcut))
valp[0].data=valp[0].data-valp[0].data[0]
valp[0].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=tcut))
tcut=4940
corr[0].trim(starttime=time_epi-timedelta(seconds=360),endtime=time_epi+timedelta(seconds=tcut))
corr[0].data=corr[0].data-corr[0].data[0]
corr[0].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=tcut))
tcut=80*60
talc[0].trim(starttime=time_epi-timedelta(seconds=360),endtime=time_epi+timedelta(seconds=tcut))
talc[0].data=talc[0].data-talc[0].data[0]
talc[0].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=tcut))



valp.write(u'/Users/dmelgar/Maule2010/wave_gauges/valp.tsun',format='SAC')
talc.write(u'/Users/dmelgar/Maule2010/wave_gauges/talc.tsun',format='SAC')
corr.write(u'/Users/dmelgar/Maule2010/wave_gauges/corr.tsun',format='SAC')


