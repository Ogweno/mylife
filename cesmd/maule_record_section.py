from numpy import genfromtxt,where
from glob import glob
from obspy.core import UTCDateTime
from obspy import read
from matplotlib import pyplot as plt
from pyproj import Geod
import matplotlib.patches as patches

rootpath=u'/Users/dmelgar/Maule2010/strong_motion/cesmd/proc/'
time_epi=UTCDateTime('2010-02-2T06:34:11Z')
tcut=180
lw=0.2
tprevious=0
sta_all=genfromtxt('/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=0,dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=[1,2])
dist=genfromtxt('/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=3)
epicenter=[-72.733002,-35.909000,35]
scale=0.1#0.05#0.18
yl=[-500,500]
xl=[0,180]
fc=[1./50,0.5]
g=Geod(ellps='WGS84')

#Rjb contour

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

for k in range(len(sta_all)):
    
    station=sta_all[k]
    print station
    
    lon=lonlat[k,0]
    lat=lonlat[k,1]
    
    #North or south
    if lat>epicenter[1]:
        n_or_s=1.0
    else:
        n_or_s=-1.0
        
    n=read(rootpath+station+'.HNN.sac')
    e=read(rootpath+station+'.HNE.sac')
    z=read(rootpath+station+'.HNZ.sac')
    
    
    #trim
    n[0].trim(starttime=time_epi,endtime=time_epi+tcut,pad=True,fill_value=0)
    e[0].trim(starttime=time_epi,endtime=time_epi+tcut,pad=True,fill_value=0)
    z[0].trim(starttime=time_epi,endtime=time_epi+tcut,pad=True,fill_value=0)
    
    
    #Scale data
    axn.plot(n[0].times(),(n[0].data/scale)+(n_or_s*dist[k]),'k',lw=lw,c='#DC143C')
    axn.set_ylim(yl)
    axn.set_xlim(xl)
    axn.set_xlabel('Seconds since OT')
    axn.set_title('North')
    axn.set_ylabel('Epicentral distance (km)')

    axe.plot(e[0].times(),(e[0].data/scale)+(n_or_s*dist[k]),'k',lw=lw,c='#6A5ACD')
    axe.set_ylim(yl)
    axe.set_xlim(xl)
    axe.set_title('East')
    axe.set_xlabel('Seconds since OT')

    axu.plot(z[0].times(),(z[0].data/scale)+(n_or_s*dist[k]),'k',lw=lw,c='#228B22')
    axu.set_ylim(yl)
    axu.set_xlim(xl)
    axu.set_xlabel('Seconds since OT')
    axu.set_title('Up')

axn.plot([30,30],yl,'k--',lw=0.8)
axn.plot([60,60],yl,'k--',lw=0.8)
axe.plot([30,30],yl,'k--',lw=0.8)
axe.plot([60,60],yl,'k--',lw=0.8)
axu.plot([30,30],yl,'k--',lw=0.8)
axu.plot([60,60],yl,'k--',lw=0.8)

#center=(5,-17)   
#axn.add_patch(patches.Rectangle(center,40,35,facecolor='w',zorder=9997))
center=(170,0)
scaleg=(9.81/2.)/scale
x=[center[0],center[0]]
y=[center[1]-scaleg/2.,center[1]+scaleg/2.]
axn.plot(x,y,lw=4,c='k',zorder=9998)
axn.annotate('50%g',xy=(135,-5),zorder=9999,fontsize=16)
                  
plt.show()