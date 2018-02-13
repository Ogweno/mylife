from matplotlib import pyplot as plt
from obspy import read
from numpy import genfromtxt
from obspy.core import UTCDateTime

#Plot acceleration
stations=genfromtxt('/Users/dmelgar/Slip_inv/wenc_2008/data/station_info/wenc.gflist',usecols=0,dtype='S')
stations=stations[5:]
path='/Users/dmelgar/Slip_inv/wenc_2008/'
time_epi=UTCDateTime('2008-05-12T06:28:04')

for k in range(len(stations)):
    sta=stations[k]
    print sta
    n=read(path+'data/waveforms/'+sta+'.vel.n')
    e=read(path+'data/waveforms/'+sta+'.vel.e')
    u=read(path+'data/waveforms/'+sta+'.vel.u')
    n.trim(starttime=time_epi)
    e.trim(starttime=time_epi)
    u.trim(starttime=time_epi)
    
    plt.figure()
    plt.subplot(311)
    plt.title(sta)
    plt.plot(n[0].times(),n[0].data)
    plt.ylabel('North (m/s)')
    plt.subplot(312)
    plt.plot(e[0].times(),e[0].data)
    plt.ylabel('East (m/s)')
    plt.subplot(313)
    plt.plot(u[0].times(),u[0].data)
    plt.ylabel('Up (m/s)')
    plt.savefig(path+'plots/'+sta+'.pdf')
 
# Plot GPS   
stations=genfromtxt('/Users/dmelgar/Slip_inv/wenc_2008/data/station_info/wenc.gflist',usecols=0,dtype='S')
stations=stations[0:5]

for k in range(len(stations)):
    sta=stations[k]
    print sta
    n=read(path+'data/waveforms/TRAKWC.NEU.'+sta+'.LC.kdisp.n')
    e=read(path+'data/waveforms/TRAKWC.NEU.'+sta+'.LC.kdisp.e')
    u=read(path+'data/waveforms/TRAKWC.NEU.'+sta+'.LC.kdisp.u')
    n.trim(starttime=time_epi)
    e.trim(starttime=time_epi)
    u.trim(starttime=time_epi)
    
    plt.figure()
    plt.subplot(311)
    plt.title(sta)
    plt.plot(n[0].times(),n[0].data)
    plt.ylabel('North (m)')
    plt.subplot(312)
    plt.plot(e[0].times(),e[0].data)
    plt.ylabel('East (m)')
    plt.subplot(313)
    plt.plot(u[0].times(),u[0].data)
    plt.ylabel('Up (m)')
    plt.savefig(path+'plots/'+sta+'.pdf')