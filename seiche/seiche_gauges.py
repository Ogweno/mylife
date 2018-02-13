from pandas import read_table, DataFrame
import datetime
import utide
from numpy import genfromtxt,zeros,mean,arange,where
from obspy import UTCDateTime
from matplotlib import pyplot as plt

f='/Users/dmelgar/seiche/alameda.csv'
lat=37+46.3/60
t0=UTCDateTime('2004-12-26')

#read pre-event
dates=genfromtxt(f,usecols=0,dtype='S')
times=genfromtxt(f,usecols=1,dtype='S')
data=genfromtxt(f,usecols=2)




#Convert to days since t0
td=zeros(len(dates))
for k in range(len(dates)):
    td[k]=(UTCDateTime(dates[k]+'T'+times[k])-t0)/86400.
    

#only include times before the eq
i=where(td<0)[0]

coef = utide.solve(td[i], data[i],
             lat=lat,
             nodal=True,
             trend=False,
             method='robust',
             conf_int='linear')
             
tide = utide.reconstruct(td, coef)




fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True, figsize=(17, 5))
c1='#0080FF'
c2='#FF4500'
c3='#3CB371'

ax0.plot(td, data, lw=0.5,c=c1)
ax0.set_ylabel('Observed (m)')
ax0.grid()

ax1.plot(td, tide['h'], lw=0.5, c=c2)
ax1.set_ylabel('Modeled (m)')
ax1.grid()

ax2.plot(td, data-tide['h'], lw=0.5,c=c3)
ax2.set_ylabel('Residual (m)')
ax2.set_xlabel('Days after September 8th, 2017')
ax2.grid()







plt.show()