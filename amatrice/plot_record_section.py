from numpy import genfromtxt,where
from glob import glob
from obspy.core import UTCDateTime
from obspy import read
from matplotlib import pyplot as plt
from pyproj import Geod
import matplotlib.patches as patches

rootpath='/Users/dmelgar/Amatrice2016/strong_motion/'
time_epi=UTCDateTime('2016-08-24T01:36:32')
tcut=30
tprevious=0
sta_all=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[1,2])
epicenter=[13.234, 42.698]
scale=0.005#0.05#0.18
yl=[-100,100]
fc=[1./50,0.5]
g=Geod(ellps='WGS84')

#Rjb contour
rjb_poly=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/Rjb.txt')

fig, axarr = plt.subplots(1, 3,figsize=(20,12))
axn=axarr[0]
axe=axarr[1]
axu=axarr[2]


def bandpass_filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt


files=glob(rootpath+'sac/*')
for k in range(len(files)):
    
    station=files[k].split('/')[-1].split('.')[0]
    print station
    
    i=where(sta_all==station)[0]
    lon=lonlat[i,0]
    lat=lonlat[i,1]
    
    #Distance to Rjb poly
    az,baz,dist=g.inv(lon,lat,epicenter[0],epicenter[1])
    dist=dist/1000.
    
    #North or south
    if lat>epicenter[1]:
        n_or_s=1.0
    else:
        n_or_s=-1.0
    
    st=read(files[k])
    accel=st[0].data
    accel=accel-accel.mean()
    if fc !=None:
        accel=bandpass_filter(accel,fc,1./st[0].stats.delta,2)
    st[0].data=accel
    
    #trim
    st[0].trim(starttime=time_epi,endtime=time_epi+tcut)
    
    chan=st[0].stats.channel
    if chan=='HLN':
        chan='HNN'
    if chan=='HLE':
        chan='HNE'
    if chan=='HLZ':
        chan='HNZ'
    
    print '...'+str(dist*n_or_s)
    
    #Scale data
    if chan=='HNN':
        axn.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
        axn.set_ylim(yl)
        axn.set_xlabel('Seconds since OT')
        axn.set_title('North')
        axn.set_ylabel('Epicentral distance (km)')
    elif chan=='HNE':
        axe.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
        axe.set_ylim(yl)
        axe.set_title('East')
        axe.set_xlabel('Seconds since OT')
    elif chan=='HNZ':
        axu.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
        axu.set_ylim(yl)
        axu.set_xlabel('Seconds since OT')
        axu.set_title('Up')
    else:
        print 'Weird channel at station '+st[0].stats.station+', '+st[0].stats.channel

center=(2,75)   
axn.add_patch(patches.Rectangle(center,5,20,facecolor='w',zorder=9997))
center=(3,85)
scaleg=(9.81/5.)/scale
x=[center[0],center[0]]
y=[center[1]-scaleg/2.,center[1]+scaleg/2.]
axn.plot(x,y,lw=4,c='k',zorder=9998)
axn.annotate('20%g',xy=(3.5,85),zorder=9999)
                  
plt.show()