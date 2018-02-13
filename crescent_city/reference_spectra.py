from numpy import savetxt,r_,where,mean,arange,zeros,array,c_
from obspy import read
from matplotlib import pyplot as plt
import nitime.algorithms as tsa
from mtspec import mtspec
from scipy.interpolate import interp1d
from run_filt import RunningMedian as rm



st=read('/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_2017.sac')
#fout='/Users/dmelgar/tidegauge_noise/crescent_city/_spectra/cres_2017_psd.txt'
st[0].data=st[0].data-mean(st[0].data)
Tw=5
Ntapers=8
psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)




#resample spectra to regular periods
period=(1./f)/60
period_interp=arange(12,180,0.001)
f=interp1d(period,psd)
psd_interp=f(period_interp)

#Smooth over
psd_median=rm(psd_interp,700)
fill1=(len(psd_interp)-len(psd_median))/2
fill2=fill1+1
psd_median=r_[zeros(fill1),array(psd_median),zeros(fill2)]


plt.figure()
plt.plot(period,psd)
plt.plot(period_interp,psd_median)
plt.xlim([0,180])
plt.ylim([0,12])
plt.xlabel('Period (min)')
plt.ylabel('PSD')
plt.show()


savetxt(fout,c_[period_interp,psd_median],fmt='%.4f',header='period(s),psd')