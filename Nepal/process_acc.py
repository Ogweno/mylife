# -*- coding: utf-8 -*-
from obspy import read
from obspy.core import UTCDateTime
from datetime import timedelta
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees,gps2DistAzimuth
from numpy import array,float64,mean
from scipy.signal import filtfilt,butter
from scipy.integrate import cumtrapz


time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=array([84.708,28.147,15]) 
delta=timedelta(seconds=0)
lon=85.336
lat=27.738
correction=38
fcorner=1./20

def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout

def hfilter(tr,fcorner,order):
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y

n=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.HNN.NQ.01')
e=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.HNE.NQ.01')
u=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.HNZ.NQ.01')

#Get p-time to site
delta=locations2degrees(lat,lon,epicenter[1],epicenter[0])
tt=getTravelTimes(delta,epicenter[2])
tp=timedelta(seconds=float64(tt[0]['time']))

#Apply shift
n[0].stats.starttime=n[0].stats.starttime-timedelta(seconds=correction)
e[0].stats.starttime=e[0].stats.starttime-timedelta(seconds=correction)
u[0].stats.starttime=u[0].stats.starttime-timedelta(seconds=correction)

#Trim
n.trim(starttime=time_epi+tp-timedelta(seconds=10),endtime=time_epi+tp+timedelta(seconds=120))
e.trim(starttime=time_epi+tp-timedelta(seconds=10),endtime=time_epi+tp+timedelta(seconds=120))
u.trim(starttime=time_epi+tp-timedelta(seconds=10),endtime=time_epi+tp+timedelta(seconds=120))

#Remove pre-event
n[0].data=n[0].data-mean(n[0].data[0:800])
e[0].data=e[0].data-mean(e[0].data[0:800])
u[0].data=u[0].data-mean(u[0].data[0:800])

#Apply gain
n[0].data=n[0].data/290000
e[0].data=e[0].data/290000
u[0].data=u[0].data/290000

##Integrate
#n[0].data=cumtrapz(n[0].data,n[0].times(),initial=0)
#e[0].data=cumtrapz(e[0].data,e[0].times(),initial=0)
#u[0].data=cumtrapz(u[0].data,u[0].times(),initial=0)
#
##Decimate
#n=stdecimate(n,5,13) ; n=stdecimate(n,4,13) ; n=stdecimate(n,2,13)
#e=stdecimate(e,5,13) ; e=stdecimate(e,4,13) ; e=stdecimate(e,2,13)
#u=stdecimate(u,5,13) ; u=stdecimate(u,4,13) ; u=stdecimate(u,2,13)
#
##High pass filter
#fnyquist=1./(2*e[0].stats.delta)
#Fc=fcorner/fnyquist
#e[0].data=hfilter(e[0].data,Fc,3)
#n[0].data=hfilter(n[0].data,Fc,3)
#u[0].data=hfilter(u[0].data,Fc,3)

#Save
#n.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.n',format='SAC')
#e.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.e',format='SAC')
#u.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.u',format='SAC')
#
##Location
#n=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.n')
#e=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.e')
#u=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.u')
#n[0].stats['sac']['stlo']=lon
#e[0].stats['sac']['stlo']=lon
#u[0].stats['sac']['stlo']=lon
#n[0].stats['sac']['stla']=lat
#e[0].stats['sac']['stla']=lat
#u[0].stats['sac']['stla']=lat

#Save
n.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.n',format='SAC')
e.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.e',format='SAC')
u.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.u',format='SAC')

n.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.n',format='SAC')
e.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.e',format='SAC')
u.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.u',format='SAC')

#Location
n=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.n')
e=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.e')
u=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.u')
n[0].stats['sac']['stlo']=lon
e[0].stats['sac']['stlo']=lon
u[0].stats['sac']['stlo']=lon
n[0].stats['sac']['stla']=lat
e[0].stats['sac']['stla']=lat
u[0].stats['sac']['stla']=lat

#Save
n.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.n',format='SAC')
e.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.e',format='SAC')
u.write(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.u',format='SAC')