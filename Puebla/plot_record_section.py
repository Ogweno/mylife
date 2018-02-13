from numpy import genfromtxt,where
from glob import glob
from obspy.core import UTCDateTime
from obspy import read
from matplotlib import pyplot as plt
from pyproj import Geod
import matplotlib.patches as patches

rootpath='/Users/dmelgar/Puebla2017/strong_motion/sac/'
time_epi=UTCDateTime('2017-09-19T18:14:37.06')
tcut=90
tprevious=0
sta_all=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/really_all_stations.txt',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/really_all_stations.txt',usecols=[1,2])
epicenter=[-98.6878,18.3044,57.52]
scale=0.3#0.05#0.18
max_dist=200
yl=[-max_dist,max_dist]
yl=[30,max_dist]
fc=[1./50,0.5]
g=Geod(ellps='WGS84')

##Rjb contour
#rjb_poly=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/Rjb.txt')

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


files=glob(rootpath+'*.sac')
for k in range(len(files)):
    
    station=files[k].split('/')[-1].split('.')[0]
    print station
    
    i=where(sta_all==station)[0]
    lon=lonlat[i,0]
    lat=lonlat[i,1]
    
    #Distance to Rjb poly
    az,baz,dist=g.inv(lon,lat,epicenter[0],epicenter[1])
    dist=dist/1000.
    
    if dist<max_dist:
        #North or south
        n_or_s=1.0
        #if lat>epicenter[1]:
        #    n_or_s=1.0
        #else:
        #    n_or_s=-1.0
        
        st=read(files[k])
        accel=st[0].data
        accel=accel-accel.mean()
        #if fc !=None:
        #    accel=bandpass_filter(accel,fc,1./st[0].stats.delta,2)
        st[0].data=accel
        
        #trim
        st[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        
        chan=files[k].split('/')[-1].split('.')[1]
        if chan=='HLN' or chan =='HNN':
            chan='HNN'
        if chan=='HLE' or chan=='HNE':
            chan='HNE'
        if chan=='HLZ' or chan=='HNZ':
            chan='HNZ'
        
        print '... '+str(dist*n_or_s)
        
        #Scale data
        if chan=='HNN':
            axn.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
            axn.set_ylim(yl)
            axn.set_xlabel('Seconds since OT',fontsize=14)
            axn.set_title('North',fontsize=14)
            axn.set_ylabel('Epicentral distance (km)',fontsize=14)
        elif chan=='HNE':
            axe.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
            axe.set_ylim(yl)
            axe.set_title('East',fontsize=14)
            axe.set_xlabel('Seconds since OT',fontsize=14)
        elif chan=='HNZ':
            axu.plot(st[0].times(),(st[0].data/scale)+(n_or_s*dist),'k',lw=0.5)
            axu.set_ylim(yl)
            axu.set_xlabel('Seconds since OT',fontsize=14)
            axu.set_title('Up',fontsize=14)
        else:
            print 'Weird channel at station '+files[k]+' :: '+st[0].stats.station+', '+st[0].stats.channel

center=(65,50)   
axn.add_patch(patches.Rectangle(center,5,20,facecolor='w',zorder=9997))
center=(65,50)
scaleg=(9.81/5.)/scale
x=[center[0],center[0]]
y=[center[1]-scaleg/2.,center[1]+scaleg/2.]
axn.plot(x,y,lw=4,c='k',zorder=9998)
axn.annotate('20%g',xy=(68,48),zorder=9999,fontsize=14)
                  
plt.show()