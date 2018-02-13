from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.event import read_events,Catalog
from obspy import read
from matplotlib import pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 18})

#Path to write waveforms to:
out='/Users/dmelgar/ElMayor2010/seed/'
t1=UTCDateTime('2010-04-04T22:40:42')
dt=120.

#phases
stime=UTCDateTime('2010-04-04T22:41:02.39')-t1
ptime=UTCDateTime('2010-04-04T22:40:53.84')-t1

#Init the client
client = Client('SCEDC')

accel_gain=262704.
vel_gain=6.273e8

st = client.get_waveforms('CI','WES','','HN_',t1,t1+dt)
a=st[1].copy()
a.data=a.data/accel_gain

st = client.get_waveforms('CI','WES','','HH_',t1,t1+dt)
v=st[2].copy()
v.data=v.data/vel_gain

d=read(u'/Users/dmelgar/ElMayor2010/proc/p494.LXN.sac')
d=d[0].copy()

#trim all to same length
d.trim(starttime=t1,endtime=t1+dt)

pltdt=90

plt.figure(figsize=(10,12))

ax=plt.subplot(311)
plt.plot(v.times(),v.data*100.,'k',lw=1)
plt.scatter(ptime,0,marker='|',c='r',s=2000,lw=3)
plt.scatter(stime,0,marker='|',c='b',s=2000,lw=3)
plt.xlim([0,pltdt])
plt.ylim([-2,2])
plt.ylabel(r'$v_z (cm/s)$')
ax.set_xticklabels([])
plt.annotate('Broadband',xy=(2,1.5))

ax=plt.subplot(312)
plt.plot(a.times(),a.data,'k')
plt.scatter(ptime,0,marker='|',c='r',s=2000,lw=3)
plt.scatter(stime,0,marker='|',c='b',s=2000,lw=3)
plt.xlim([0,pltdt])
plt.ylabel(r'$a_n (m/s^2)$')
plt.ylim([-2,3])
ax.set_xticklabels([])
plt.annotate('Strong motion',xy=(2,2.5))

ax=plt.subplot(313)
plt.plot(d.times(),d.data,'k')
plt.scatter(ptime,0,marker='|',c='r',s=2000,lw=3)
plt.scatter(stime,0,marker='|',c='b',s=2000,lw=3)
plt.xlim([0,pltdt])
plt.ylim([-0.7,0.3])
plt.ylabel(r'$d_n (m)$')
plt.xlabel('Seconds after OT')
plt.annotate('GPS',xy=(2,0.16))

plt.subplots_adjust(left=0.1,right=0.98,bottom=0.1,top=0.98,hspace=0.08)