from mtspec import mtspec
from numpy import diff,where,c_,savetxt
from obspy import read
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

stations=['M18B','M10B','G34B','M13B','FS20B','FS13B']
path='/Users/dmelgar/tidegauge_noise/APG/'

#latitudes=[45.853,42.682,39.333,48.796]

Tw=5
Ntapers=8

plt.figure()

spike_max=20
vmax=0.3


for ksta in range(len(stations)):

    sta=stations[ksta]

    st=read(path+sta+'_2013_05.sac')
    raw_data=st[0].data-st[0].data[0]
    
    #The time I always want to resample to
    tabs=st[0].times()
    
    #remove big obvious peaks
    i=where(abs(raw_data)<spike_max)[0]
    t_no_big=tabs[i]
    d_no_big=raw_data[i]
    f=interp1d(t_no_big,d_no_big)
    dclean=f(tabs)
    
    #Look for spikes (more like weird boxcars)
    v=diff(dclean)
    i=where(abs(v)>vmax)[0]
    dclean2=dclean.copy()
    for k in range(len(i)/2):
        eta1=dclean[i[2*k]]
        eta2=dclean[i[2*k+1]+1]
        t1=tabs[i[2*k]]
        t2=tabs[i[2*k+1]+1]
        slope=(eta2-eta1)/(t2-t1)
        inter=eta1-slope*t1
        dclean2[i[2*k]:i[2*k+1]+1]=slope*tabs[i[2*k]:i[2*k+1]+1]+inter
        
    #finally remove mean
    depth=dclean2.mean()
    st[0].data=dclean2-depth
    
    #time in days
    td=st[0].times()/86400.
    
    psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)
    
    period=(1./f)/60
    
    plt.loglog(period,psd,label=sta)
    
    out=c_[period[1:],psd[1:]]
    fout='/Users/dmelgar/tidegauge_noise/APG/spectra/'+sta+'_2017_05.spec'
    savetxt(fout,out,fmt='%.5e',header='period (min), psd')
    
plt.legend()
plt.grid()
plt.xlabel('Period')
plt.ylabel('PSD')
plt.xlim([1,1600])
plt.show() 
    