from mtspec import mtspec
from obspy import read
from matplotlib import pyplot as plt
import utide
from numpy import c_,savetxt

stations=['46404','46407','46411','46419']
latitudes=[45.853,42.682,39.333,48.796]
months=['06','06','06','07']
Tw=5
Ntapers=8

plt.figure()

for k in range(len(stations)):
    sta=stations[k]
    lat=latitudes[k]
    mo=months[k]
    st=read('/Users/dmelgar/tidegauge_noise/DART/sac/'+sta+'_2017_'+mo+'.sac')
    depth=st[0].data.mean()
    st[0].data=st[0].data-depth
    
    #time in days
    td=st[0].times()/86400.
    
    #get tides
    #coef = utide.solve(td, st[0].data,lat=lat,nodal=True,trend=False,method='robust',conf_int='linear')
    #tide = utide.reconstruct(td, coef)
    #tide=tide['h']
    ##data=st[0].data-tide
    data=st[0].data
    
    
    psd, f, jackknife, _, _ = mtspec(data=data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)
    
    period=(1./f)/60
    
    plt.loglog(period,psd,label=sta+'('+str(int(depth))+'m)')
    
    out=c_[period[1:],psd[1:]]
    fout='/Users/dmelgar/tidegauge_noise/DART/spectra/'+sta+'_2017_'+mo+'.spec'
    savetxt(fout,out,fmt='%.5e',header='period (min), psd')
    
plt.legend()
plt.grid()
plt.xlabel('Period')
plt.ylabel('PSD')
plt.xlim([30,1600])
plt.show() 
    