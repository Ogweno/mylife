from pandas import read_table, DataFrame
import datetime
import utide
from numpy import genfromtxt,zeros,mean,arange
from obspy import UTCDateTime
from matplotlib import pyplot as plt

f='/Users/dmelgar/Chiapas2017/post_seismic/tide_gauges/sali_pre_earthquake.txt'
lat=16.168436
t0=UTCDateTime('2017-09-08')

#post earthquake data
f2='/Users/dmelgar/Chiapas2017/post_seismic/tide_gauges/sali_after_earthquake.txt'
delta=60/86400.
days_pos=90

#read pre-event
dates=genfromtxt(f,usecols=0,dtype='S')
times=genfromtxt(f,usecols=1,dtype='S')
data=genfromtxt(f,usecols=2)


#read post-event
dates_post=genfromtxt(f2,usecols=0,dtype='S')
times_post=genfromtxt(f2,usecols=1,dtype='S')
data_post=genfromtxt(f2,usecols=2)



#Convert to days since t0
td=zeros(len(dates))
for k in range(len(dates)):
    td[k]=(UTCDateTime(dates[k]+'T'+times[k])-t0)/86400.
    
#days after t0 for post
td_post=zeros(len(dates_post))
for k in range(len(dates_post)):
    td_post[k]=(UTCDateTime(dates_post[k]+'T'+times_post[k])-t0)/86400.
    
#resample to regular interval
    
    
coef = utide.solve(td, data,
             lat=lat,
             nodal=True,
             trend=False,
             method='robust',
             conf_int='linear')
             
tide = utide.reconstruct(td, coef)

#Predcit for at the post event times
tide_predict=utide.reconstruct(td_post, coef)




fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True, figsize=(17, 5))
c1='#0080FF'
c2='#FF8C00'

ax0.plot(td, data, c=c1)
ax0.plot(td_post, data_post, c=c2)

ax1.plot(td, tide['h'], alpha=0.5, c=c1)
ax1.plot(td_post, tide_predict['h'], alpha=0.5, c=c2)

ax2.plot(td, data-tide['h'], alpha=0.5,c=c1)
ax2.plot(td_post, data_post-tide_predict['h'], alpha=0.5,c=c2)






plt.show()