from numpy import genfromtxt,arange,where
from scipy.interpolate import interp1d
from mtspec import mtspec
from matplotlib import pyplot as plt

g=genfromtxt('/Users/dmelgar/Tsunamis/crescent_city/_output/gauge00101.txt')
#g=genfromtxt('/Users/dmelgar/Tsunamis/crescent_city/20km_1e-2/_output/gauge00001.txt')
#g=genfromtxt('/Users/dmelgar/Tsunamis/crescent_city/20km_48hrs_1e-3/_output/gauge00001.txt')
#g=genfromtxt('/Users/dmelgar/Tsunamis/crescent_city/20km_48hrs_1e-3_south/_output/gauge00101.txt')

i=where(g[:,5]>1)[0]
g[i,5]=0
plt.figure()
plt.plot(g[:,1],g[:,5])
plt.show()

dt=10.
Tw=3
Ntapers=8
tg=g[:,1]
eta_unfilt=g[:,5]

#Get refernce spectra
ref_psd=genfromtxt(u'/Users/dmelgar/tidegauge_noise/data/cres/spectra/cres_2017_psd.txt')


#filter
fcorner=[1./(200*60),1./(1*60)]
def bandpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt



#resample tor egular itnerval
t=arange(0,tg.max(),dt)
f=interp1d(tg,eta_unfilt)
eta_unfilt=f(t)

#filter
eta=eta_unfilt#bandpass(eta_unfilt,fcorner,1./dt,2)

#Get psd
psd, f = mtspec(data=eta, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=len(eta), statistics=False)
psd_u, f = mtspec(data=eta_unfilt, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=len(eta), statistics=False)
period=1./f/60


plt.figure()
plt.plot(period,psd/psd.max())
#plt.plot(period,psd)
plt.plot(ref_psd[:,0],ref_psd[:,1]/ref_psd[:,1].max())
#plt.plot(ref_psd[:,0],ref_psd[:,1])
plt.xlim([12,100])
plt.xlabel('Period (min)')
plt.legend(['Gaussian source','Reference PSD'])
plt.ylabel('PSD')
plt.grid()
plt.show()