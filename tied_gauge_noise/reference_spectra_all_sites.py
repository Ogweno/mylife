from numpy import savetxt,r_,where,mean,arange,zeros,array,c_,genfromtxt,median,mean
from obspy import read
from matplotlib import pyplot as plt
from mtspec import mtspec
from scipy.interpolate import interp1d
from run_filt import RunningMedian as rm



stations=genfromtxt('/Users/dmelgar/tidegauge_noise/station_info/stations.txt',usecols=1,dtype='S')
Tw=25
Ntapers=15
max_period=120
Nchunks=24

for k in range(len(stations)):
    
    sta=stations[k]
    print sta
    
    
    st=read('/Users/dmelgar/tidegauge_noise/data/'+sta+'/sac/'+sta+'_2017.sac')
    fout='/Users/dmelgar/tidegauge_noise/data/'+sta+'/spectra/'+sta+'_2017_psd.txt'
    st[0].data=st[0].data-mean(st[0].data)
    
    #Get the spectra in chunks
    tmax=st[0].times()[-1]
    N=st[0].stats.npts
    segments=arange(0,N,int(N/Nchunks))

    plt.figure(figsize=(12,3.8))
    for kspec in range(Nchunks):
        N1=segments[kspec]
        if kspec<Nchunks-1:
            N2=segments[kspec+1]
        else:
            N2=N
    
        psd, f, jackknife, _, _ = mtspec(data=st[0].data[N1:N2], delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)
        if kspec==0:
            period=(1./f)/60
            i=where((period>12) & (period<max_period))[0]
            period=period[i] #reference period
            psd_out=psd[i]
            plt.plot(period,psd_out)
            
        if kspec>0:
            #resample to reference period
            current_period=(1./f[1:])/60
            current_psd=psd[1:]
            f=interp1d(current_period,current_psd)
            psd_interp=f(period)
            psd_out=c_[psd_out,psd_interp]
            plt.plot(period,psd_interp)
        
    psd_out=median(psd_out,axis=1)
    plt.title(sta)
    plt.plot(period,psd_out,c='k',lw=1.5)
    plt.xlim([12,period.max()])
    plt.ylim([0,psd_out.max()*1.2])
    plt.xlabel('Period (min)')
    plt.ylabel('PSD')
    plt.grid()

    plt.savefig('/Users/dmelgar/tidegauge_noise/plots/spectra/'+sta+'.png')
    plt.close()

    savetxt(fout,c_[period,psd_out],fmt='%.4f',header='period(s),psd')