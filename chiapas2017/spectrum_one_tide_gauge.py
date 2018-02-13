import matplotlib.pyplot as plt
from obspy import read
from numpy import genfromtxt,arange,where,amax
from scipy.interpolate import interp1d
from obspy import UTCDateTime,Stream,Trace
from mtspec import mtspec
from datetime import timedelta
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

rcParams.update({'font.size': 16})
dt=60.
time_epi=UTCDateTime('2017-09-08T04:49:21')
trim_time=timedelta(seconds=48*3600)
stations=['ptan','huat','sali','pchi']
gauges=['01001','01002','01003','01004']
ylims_obs=[0.3,1.0,1.4,2]
ylims_syn=[0.3,1.0,1.4,1.2]
flims=[0.3,10]
#tmax=[24,24,24.,24.]
tmax=[48,48,48.,48.]
plims=[[0.0005,5],[0.002,90],[0.02,600],[0.0002,3000]]
Tband=[12,4,10,4]

path=u'/Users/dmelgar/Chiapas2017/plots/spec/'


f=plt.figure(figsize=(7.6,7))
obsc='#1E90FF'
sync='#DAA520'



def lowpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/fnyquist,'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt

#Corner frequency (10min period lowpass)
fc=[1./7200,10/3600.]


for k in range(len(stations)):
    obs=read(u'/Users/dmelgar/Chiapas2017/tsunami/long/'+stations[k]+'.sac')
    #gauge=genfromtxt('/Users/dmelgar/Tsunamis/tehuantepec_hf_all/_output/gauge'+gauges[k]+'.txt')
    #gauge=genfromtxt('/Users/dmelgar/Tsunamis/tehuantepec_hf_deepgauges/_output/gauge'+gauges[k]+'.txt')
    #gauge=genfromtxt('/Users/dmelgar/Tsunamis/tehuantepec_hf_deepgauges_v2/_output/gauge'+gauges[k]+'.txt')
    gauge=genfromtxt('/Users/dmelgar/Tsunamis/tehuantepec_48hrs/_output/gauge'+gauges[k]+'.txt')
    
    if k==2:
        gauge=genfromtxt('/Users/dmelgar/Tsunamis/tehuantepec_hf_all/_output/gauge'+gauges[k]+'.txt')
        #shoal=(abs(1000)/100)**0.25 
        shoal=1
    else:
        shoal=1.
    
    #Trim data to one day
    obs[0].trim(starttime=time_epi,endtime=time_epi+trim_time)
    obs[0].stats.delta=60
    
    #Resample synthetics tor egular itnerval
    tsyn=obs[0].times()[0:-1]
    i=where(tsyn<24*3600)[0]
    tsyn=tsyn[i]
    
    f=interp1d(gauge[:,1],gauge[:,5])
    syn=f(tsyn)*shoal
    
    # sample rate
    delta=obs[0].stats.delta
    
    #Fitler the thingamajigs
    obs[0].data=lowpass(obs[0].data,fc,1./delta,4)
    syn=lowpass(syn,fc,1./delta,4)
    
    #get spectra
    i=where(obs[0].times()<tmax[k]*3600)[0]
    Sobs, fobs = mtspec(
        data=obs[0].data[i], delta=delta, time_bandwidth=3.5,
        number_of_tapers=5, nfft=len(obs[0].data[i]))
    
    i=where(tsyn<tmax[k]*3600)[0]    
    Ssyn, fsyn = mtspec(
        data=syn[i], delta=delta, time_bandwidth=Tband[k],
        number_of_tapers=5, nfft=len(syn[i]), adaptive=True)
    
    #Use grid spec to define the axes    
    gs = gridspec.GridSpec(8, 2,width_ratios=[2,1])
    ax1=plt.subplot(gs[2*k,0])
    ax2=plt.subplot(gs[2*k+1,0])
    ax3=plt.subplot(gs[2*k:2*k+2,1])

    #Define the axes    
    #ax1 = plt.subplot2grid((8, 2), (2*k, 0))
    #ax2 = plt.subplot2grid((8, 2), (2*k+1, 0))
    #ax3 = plt.subplot2grid((8, 2), (2*k, 1),rowspan=2)
    
    
    
    #noise
    #pos = [0.05,0.6,0.4,0.38]
    #ax1.set_position(pos)
    ax1.plot(obs[0].times()/3600.,obs[0].data,c=obsc,lw=0.5)
    ax1.set_xlim([0,24])
    ax1.set_ylim([-ylims_obs[k],ylims_obs[k]])
    ax1.grid()
    ax1.set_ylabel(stations[k])
    ax1.yaxis.set_ticklabels([])
    
    #pos = [0.05,0.14,0.4,0.38]
    #ax2.set_position(pos)
    ax2.plot(tsyn/3600.,syn,c=sync,lw=0.5)
    ax2.set_xlim([0,24])
    ax2.set_ylim([-ylims_syn[k],ylims_syn[k]])
    ax2.grid()
    ax2.yaxis.set_ticklabels([])
    
    
    #pos = [0.76,0.3,0.22,0.48]
    #ax3.set_position(pos)
    ax3.set_xlim(flims)
    ax3.set_ylim(plims[k])
    ax3.loglog(fobs*3600,Sobs,c=obsc)
    ax3.loglog(fsyn*3600,Ssyn,c=sync)
    ax3.grid(which='both')
    ax3.yaxis.tick_right()
    
    #ax3.set_ylim([3e-6,500])
    
    if k==0:
        ax1.xaxis.set_ticklabels([])
        ax2.xaxis.set_ticklabels([])
        ax3.xaxis.set_ticklabels([])
        ax1.plot([-1,-1],[0,0],c=sync)
        ax1.legend(['obs','syn'],bbox_to_anchor=(0.1, 0.8),ncol=2,frameon=False)
        y=0.1
        ax1.annotate(str(round(amax(obs[0].data),3)),xy=(20,y),fontsize=12,color=obsc)
        ax2.annotate(str(round(amax(syn),3)),xy=(20,y),fontsize=12,color=sync)
        
    if k==1:
        ax1.xaxis.set_ticklabels([])
        ax2.xaxis.set_ticklabels([])
        ax3.xaxis.set_ticklabels([])
        y=0.3
        ax1.annotate(str(round(amax(obs[0].data),3)),xy=(20,y),fontsize=12,color=obsc)
        ax2.annotate(str(round(amax(syn),3)),xy=(20,y),fontsize=12,color=sync)
        
    if k==2:
        ax1.xaxis.set_ticklabels([])
        ax2.xaxis.set_ticklabels([])
        ax3.xaxis.set_ticklabels([])
        y=0.5
        ax1.annotate(str(round(amax(obs[0].data),3)),xy=(20,y),fontsize=12,color=obsc)
        ax2.annotate(str(round(amax(syn),3)),xy=(20,y),fontsize=12,color=sync)
    
    
    if k==3:
        ax1.xaxis.set_ticklabels([])
        ax2.set_xlabel('Hours after OT')
        ax3.set_xlabel('cph')
        y1=0.6
        y2=0.3
        ax1.annotate(str(round(amax(obs[0].data),3)),xy=(20,y1),fontsize=12,color=obsc)
        ax2.annotate(str(round(amax(syn),3)),xy=(20,y2),fontsize=12,color=sync)
        
    
plt.subplots_adjust(hspace=0,wspace=0.02)
    

plt.show()

#plt.savefig(path+station+'.png')