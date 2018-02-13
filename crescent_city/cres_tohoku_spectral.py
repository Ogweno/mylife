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

#st=read('/Users/dmelgar/tidegauge_noise/eureka/eure_2017.sac')
#vmin=-1.5
#vmax=0.5
#sealim=0.15
#specmax=2.5

st=read('/Users/dmelgar/tidegauge_noise/crescent_city/cres_03_2011.sac')
vmin=0
vmax=3
sealim=2.5
specmax=500

st[0].data=st[0].data-mean(st[0].data)
Tw=5
Ntapers=12

Ndays=0.4
Nsamples=int(Ndays*6*60)
dt=int(round(st[0].stats.delta))


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
    

#### WHole spectrum

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



#go to period
period_spec=1./fspec[1:]


#Fitler
fcorner=1./(2*3600)
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
st_filt=st.copy()
st_filt[0].data=highpass(st[0].data,fcorner,1./st[0].stats.delta,2)



plt.figure(figsize=(16,6))


#spectrogram
ax=plt.axes([0.05,0.1,0.68,0.53])
obj=ax.pcolormesh(time/86400,period_spec/60,log10(psd_out[1:,:]),vmin=vmin,vmax=vmax,cmap=plt.cm.jet)
#cb=plt.colorbar(obj)
#cb.set_label('dB')
ax.set_ylim([12,100])
ax.set_xlim([0,30])
ax.set_ylabel('Period (min)')
ax.set_xlabel('Days')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
xmajorLocator = MultipleLocator(5)
xminorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)


#Waveforms
ax=plt.axes([0.05,0.65,0.68,0.24])
ax.xaxis.tick_top()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(st_filt[0].times()/86400,st_filt[0].data,c='#303030',lw=0.5)
ax.set_xlim([0,30])
ax.set_ylim([-sealim,sealim])
ax.set_ylabel('Sea level (m)')
ax.set_xlabel('Days')
ax.xaxis.set_label_position("top")
xmajorLocator = MultipleLocator(5)
xminorLocator = MultipleLocator(1)
ymajorLocator = MultipleLocator(1)
yminorLocator = MultipleLocator(0.2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)



#Whole spectrum
ax=plt.axes([0.74,0.1,0.2,0.53])
ax.yaxis.tick_right()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(psd_interp,period_interp,lw=1,c='#303030')
ax.set_ylim([12,100])
ax.set_xlim([0,specmax])
xmajorLocator = MultipleLocator(100)
xminorLocator = MultipleLocator(20)
ymajorLocator = MultipleLocator(20)
yminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_label_position("right")
ax.set_xlabel(r'PSD ($m^2\times\; min$)')
ax.set_ylabel('Period (min)')







# Analysis of only tohoku bit and compare to whole year
t1=UTCDateTime('2011-03-11T00:00:00')
t2=UTCDateTime('2011-03-14T00:00:00')
st[0].trim(starttime=t1,endtime=t2)
Tw=5
Ntapers=8
psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)
period=1./f/60

#get the whole year
st2=read('/Users/dmelgar/tidegauge_noise/crescent_city/cres_2017.sac')
psd2, f2, jackknife, _, _ = mtspec(data=st2[0].data, delta=st2[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st2[0].stats.npts, statistics=True)
period2=1./f2/60
#resample spectra to regular periods
period_interp=arange(12,120,0.001)
finterp=interp1d(period2,psd2)
psd_interp=finterp(period_interp)

#Smooth over
psd_median=rm(psd_interp,500)
fill1=(len(psd_interp)-len(psd_median))/2
fill2=fill1+1
psd_median=r_[zeros(fill1),array(psd_median),zeros(fill2)]


plt.figure(figsize=(5,5))
ax=plt.axes([0.15,0.1,0.83,0.8])

ax.plot(period,psd,lw=1,c='#303030',label='Tohoku-oki tsunami')
ax.plot(period_interp,psd_median*100,lw=1,c='#FF4500',label='200 x Noise')
ax.legend()

ax.set_xlim([12,100])
ax.set_ylim([0,1400])
ymajorLocator = MultipleLocator(100)
yminorLocator = MultipleLocator(20)
xmajorLocator = MultipleLocator(20)
xminorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.set_ylabel(r'PSD ($m^2\times\; min$)')
ax.set_xlabel('Period (min)')





plt.show()