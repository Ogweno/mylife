from numpy import genfromtxt,r_,where,mean,log10,array
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
import nitime.algorithms as tsa
from matplotlib.ticker import MultipleLocator
from mtspec import mtspec
from numpy import arange,zeros
from scipy.interpolate import interp1d
from run_filt import RunningMedian as rm
from obspy import UTCDateTime

#st=read('/Users/dmelgar/tidegauge_noise/eureka/eure_2017.sac')
#vmin=-1.5
#vmax=0.5
#sealim=0.15
#specmax=2.5

st=read(u'/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_2017.sac')
stp=read(u'/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_pres_2017.mseed')
#stp=read(u'/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_temp_2017.mseed')
vmin=-0.5
vmax=1
vmin_pres=0.5
vmax_pres=3.5
sealim=0.25
preslim=5
specmax=11
specmax_pres=50000
Tw=5
Ntapers=8
Ndays=3
Nsamples=Ndays*6*60
dt=int(round(st[0].stats.delta))



#trim
st[0].trim(starttime=UTCDateTime('20170301T00:00:00'))
stp.trim(starttime=UTCDateTime('20170301T00:00:00'))

#demean
st[0].data=st[0].data-mean(st[0].data)
stp[0].data=stp[0].data-mean(stp[0].data)



#tide gauge spectrogram

#Make time variable
time=arange(Nsamples*dt/2,st[0].times()[-1],dt*Nsamples/2)
Nslices=len(time)
print Nslices

#get spectrogram
sample1=0
sample2=sample1+Nsamples

ktime=0
for ktime in range(0,Nslices):
    # get psd of only that slice   
    psd, fspec, jackknife, _, _ = mtspec(data=st[0].data[sample1:sample2], delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=Nsamples, statistics=True)
    
    #initalize output variables
    if ktime==0:
        psd_out=zeros((len(fspec),Nslices))
        
    #place spectra inc orrect palce
    psd_out[:,ktime]=psd
    
    #Update counters 50% overlap fixed for now
    sample1+=Nsamples-Nsamples/2
    sample2+=Nsamples-Nsamples/2

#go to period
period_spec=1./fspec[1:]    
    
 
    
       
#pressure spectrogram

#Make time variable
time_pres=arange(Nsamples*dt/2,stp[0].times()[-1],dt*Nsamples/2)
Nslices=len(time_pres)
print Nslices

#get spectrogram
sample1=0
sample2=sample1+Nsamples

ktime=0
for ktime in range(0,Nslices):
    # get psd of only that slice   
    psd_pres, fspec_pres, jackknife, _, _ = mtspec(data=stp[0].data[sample1:sample2], delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=Nsamples, statistics=True)
    
    #initalize output variables
    if ktime==0:
        psd_out_pres=zeros((len(fspec),Nslices))
        
    #place spectra inc orrect palce
    psd_out_pres[:,ktime]=psd_pres
    
    #Update counters 50% overlap fixed for now
    sample1+=Nsamples-Nsamples/2
    sample2+=Nsamples-Nsamples/2   
    
#go to period
period_spec_pres=1./fspec_pres[1:]       
    
    
    

#### WHole spectrum tide gauge

Tw=5
Ntapers=8
psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)

#resample spectra to regular periods
period=(1./f)/60
period_interp=arange(12,120,0.001)
f=interp1d(period,psd)
psd_interp=f(period_interp)

#Smooth over
psd_median=rm(psd_interp,700)
fill1=(len(psd_interp)-len(psd_median))/2
fill2=fill1+1
psd_median=r_[zeros(fill1),array(psd_median),zeros(fill2)]

##################################




#### WHole spectrum pressure

Tw=5
Ntapers=8
psd_pres, f_pres, jackknife, _, _ = mtspec(data=stp[0].data, delta=stp[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=stp[0].stats.npts, statistics=True)

#resample spectra to regular periods
period=(1./f_pres)/60
period_interp_pres=arange(12,48*60,0.01)
f=interp1d(period,psd_pres)
psd_interp_pres=f(period_interp_pres)

#Smooth over
psd_median_pres=rm(psd_interp_pres,1500)
fill1=(len(psd_interp_pres)-len(psd_median_pres))/2
fill2=fill1+1
psd_median_pres=r_[zeros(fill1),array(psd_median_pres),zeros(fill2)]

##################################





#Fitler

def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt

fcorner=1./(2*3600)
st_filt=st.copy()
st_filt[0].data=highpass(st[0].data,fcorner,1./st[0].stats.delta,2)

fcorner=1./(48*3600)
st_filt_pres=stp.copy()
st_filt_pres[0].data=highpass(stp[0].data,fcorner,1./st[0].stats.delta,2)









#### Full spectrogram figure

plt.figure(figsize=(16,8))


#spectrogram
ax=plt.axes([0.05,0.51,0.68,0.28])
obj=ax.pcolormesh(58+time/86400,period_spec/60,log10(psd_out[1:,:]),vmin=vmin,vmax=vmax,cmap=plt.cm.jet)
ax.set_ylim([12,100])
ax.set_xlim([58,364])
ax.set_ylabel('Period (min)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
xmajorLocator = MultipleLocator(50)
xminorLocator = MultipleLocator(10)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.xaxis.set_ticks([])


#Waveforms
ax=plt.axes([0.05,0.8,0.68,0.12])
ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(58+st_filt[0].times()/86400,st_filt[0].data,c='#303030',lw=0.5)
ax.set_xlim([58,364])
ax.set_ylim([-sealim,sealim])
ax.set_ylabel('Sea level (m)')
ax.set_xlabel('Days')
ax.xaxis.set_label_position("top")
xmajorLocator = MultipleLocator(50)
xminorLocator = MultipleLocator(10)
ymajorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.05)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)



#Whole spectrum
ax=plt.axes([0.74,0.51,0.2,0.28])
ax.yaxis.tick_right()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(psd_interp,period_interp,lw=0.3,c='#303030')
ax.plot(psd_median,period_interp,lw=1,c='#FF4500')
ax.set_ylim([12,100])
ax.set_xlim([0,specmax])
xmajorLocator = MultipleLocator(2)
xminorLocator = MultipleLocator(0.5)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_label_position("right")
ax.set_xlabel(r'PSD ($m^2\times\; min$)')
ax.set_ylabel('Period (min)')



#### Pressure

#spectrogram
ax=plt.axes([0.05,0.07,0.68,0.28])
obj=ax.pcolormesh(58+time/86400,period_spec_pres/360,log10(psd_out_pres[1:,:]),vmin=vmin_pres,vmax=vmax_pres,cmap=plt.cm.jet)
ax.set_ylim([0.2,24])
ax.set_xlim([58,364])
ax.set_ylabel('Period (days)')
ax.set_xlabel('Days')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
xmajorLocator = MultipleLocator(50)
xminorLocator = MultipleLocator(10)
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)


#Waveforms
ax=plt.axes([0.05,0.36,0.68,0.12])
ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(58+st_filt_pres[0].times()/86400,st_filt_pres[0].data,c='#303030',lw=0.5)
ax.set_xlim([58,364])
ax.set_ylim([-preslim,preslim])
ax.set_ylabel('Pressure (mb)')
ax.xaxis.set_label_position("top")
xmajorLocator = MultipleLocator(50)
xminorLocator = MultipleLocator(10)
ymajorLocator = MultipleLocator(0.5)
yminorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.xaxis.set_ticks([])


#Whole spectrum
ax=plt.axes([0.74,0.07,0.2,0.28])
ax.yaxis.tick_right()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(1./f_pres/(60*60),psd_pres,lw=0.5,c='#303030')
ax.set_ylim([0.2,48])
ax.set_xlim([0,specmax_pres])
xmajorLocator = MultipleLocator(5000)
xminorLocator = MultipleLocator(1000)
ymajorLocator = MultipleLocator(5)
yminorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_label_position("right")
ax.set_xlabel(r'PSD ($mb^2\times\; min$)')
ax.set_ylabel('Period (daysn)')



plt.show()