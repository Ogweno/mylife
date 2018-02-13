from matplotlib import pyplot as plt
from numpy import genfromtxt,arange,array,zeros,mean
from obspy import read,Stream
from obspy.core import UTCDateTime
from matplotlib.ticker import AutoMinorLocator
import matplotlib
from datetime import timedelta
from mudpy.forward import lowpass

matplotlib.rcParams['font.size'] = 16
from obspy.core.util.geodetics import gps2DistAzimuth


outfile='/Users/dmelgar/code/GMT/coquimbo/wiggles.txt'
path=u'/Users/dmelgar/Coquimbo2015/strong_motion/proc/'
lonlat=genfromtxt('/Users/dmelgar/Coquimbo2015/strong_motion/strong_motion_plot.sta',usecols=[0,1])
sta=genfromtxt('/Users/dmelgar/Coquimbo2015/strong_motion/strong_motion_plot.sta',usecols=[2],dtype='S')
time_epi=UTCDateTime('2015-09-16T22:54:33')
epicenter=array([-71.654,-31.570,29.8])
ampl_adjust=10
t_adjust=0.005
corner=0.2

st=Stream()
st2=Stream()
st3=Stream()
max_ampl=0
d=zeros(len(lonlat))
for k in range(len(sta)):
    st+=read(path+sta[k]+'.HNN.sac')
    st2+=read(path+sta[k]+'.HNE.sac')
    st3+=read(path+sta[k]+'.HNE.sac')
    if st[k].data.max>max_ampl:
        max_ampl=st[k].data.max()
    st[k].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=200))
    st[k].data=st[k].data-mean(st[k].data[0:50])
    st2[k].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=200))
    st2[k].data=st2[k].data-mean(st2[k].data[0:50])
    st3[k].trim(starttime=time_epi,endtime=time_epi+timedelta(seconds=200))
    st3[k].data=st3[k].data-mean(st3[k].data[0:50])
    #Get PGA
    aabs=max((st[k].data**2+st2[k].data**2+st3[k].data**2)**0.5)/9.81*100
    print sta[k]+', PGA = '+str(aabs)+'%g'
    #st[k].data=lowpass(st[k].data,corner,1./st[k].stats.delta,2)
    d[k],az,baz=gps2DistAzimuth(epicenter[1],epicenter[0],lonlat[k,1],lonlat[k,0])
    d[k]=d[k]/1000
    if lonlat[k,1]<epicenter[1]: #It south
        d[k]=-d[k]
 
#Write to file
tout=array([])
aout=array([])
f=open(outfile,'w')
for k in range(len(sta)):
    t=st[k].times()
    t=(t*t_adjust)+lonlat[k,0]
    offset=lonlat[k,1]
    a=((st[k].data/1.0)/ampl_adjust)+offset
    f.write('>\n')
    for kdat in range(len(t)):
        f.write('%.6f\t%.6f\n' % (t[kdat],a[kdat]))
f.close()
    

plt.figure()      
for k in range(len(sta)):
    t=st[k].times()
    t=(t*t_adjust)+lonlat[k,0]
    offset=lonlat[k,1]
    a=(st[k].data/abs(st[k].data.max())/ampl_adjust)+offset
    plt.plot(t,a,lw=1,c='k')
plt.show()                             
                             
#fig, axarr = plt.subplots(2,1,figsize=(9,8)) 
#for k in range(len(sta)):
#    t=st[k].times()
#    delta=st[k].stats.starttime-time_epi
#    t=t+delta
#    offset=d[k]
#    a=(st[k].data/abs(st[k].data.max())/ampl_adjust)+offset
#    if d[k]>0:
#        ax=axarr[0]
#    else:
#        ax=axarr[1]
#    ax.plot(t,a,lw=1,c='k')
#axarr[0].set_xlim([0,200])
#axarr[1].set_xlim([0,200])
#
##plt.ylabel('Latitude')
##plt.xlabel('Seconds after OT')
#minorLocator = AutoMinorLocator()
#ax.xaxis.set_minor_locator(minorLocator)
#ax.tick_params(which='both', width=1)
#ax.tick_params(which='major', length=8)
#ax.tick_params(which='minor', length=3)    
#plt.show()