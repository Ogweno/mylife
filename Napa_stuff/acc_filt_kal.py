from obspy import read
from numpy import mean,log10,r_,diff
from scipy.integrate import cumtrapz
from scipy.signal import butter,filtfilt
from matplotlib import pyplot as plt
import nitime.algorithms as tsa
from datetime import timedelta

filt=False
low_corner=1.0

ae=read(u'/Users/dmelgar/Napa2014/acc/_trim/NC.N014.HNE.sac')
e=read(u'/Users/dmelgar/Napa2014/GPS/trim/p198.LXE.sac')
dke=read(u'/Users/dmelgar/Napa2014/kal/NC.N014.HXE.sac')
vke=read(u'/Users/dmelgar/Napa2014/kal/NC.N014.HYE.sac')

#Corrections to a
ae[0].data=ae[0].data-mean(ae[0].data[0:7000])
ae.trim(starttime=ae[0].stats.starttime+timedelta(seconds=38))
#Corrections to e
e.trim(starttime=ae[0].stats.starttime,endtime=ae[0].stats.endtime)
e[0].data=e[0].data-mean(e[0].data[167:205])
e.trim(starttime=e[0].stats.starttime+timedelta(seconds=38))
#Corrections to kal
dke.trim(starttime=dke[0].stats.starttime+timedelta(seconds=38))
vke.trim(starttime=vke[0].stats.starttime+timedelta(seconds=38))

#Apply filter
if filt==True:
    fnyquist=1./(2*ae[0].stats.delta)
    b, a = butter(4, low_corner/(fnyquist),'high')
    ae[0].data=filtfilt(b,a,ae[0].data)
else:
    low_corner=0

#Integrate 
vae=ae.copy()
dae=ae.copy()
vae[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
dae[0].data=cumtrapz(vae[0].data,vae[0].times(),initial=0)
#Differentiate
ve=e.copy()
ve[0].data=r_[0,diff(e[0].data)]/ve[0].stats.delta


#Spectrums of things
fdae, Pdae, nu = tsa.multi_taper_psd(dae[0].data,Fs=1./dae[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
fdke, Pdke, nu = tsa.multi_taper_psd(dke[0].data,Fs=1./dke[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
fe, Pe, nu = tsa.multi_taper_psd(e[0].data,Fs=1./e[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
Pdae=10*log10(Pdae)
Pdke=10*log10(Pdke)
Pe=10*log10(Pe)
fvae, Pvae, nu = tsa.multi_taper_psd(vae[0].data,Fs=1./vae[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
fvke, Pvke, nu = tsa.multi_taper_psd(vke[0].data,Fs=1./vke[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
fve, Pve, nu = tsa.multi_taper_psd(ve[0].data,Fs=1./ve[0].stats.delta,adaptive=False,jackknife=False,low_bias=False)
Pvae=10*log10(Pvae)
Pvke=10*log10(Pvke)
Pve=10*log10(Pve)

#plot
plt.figure()
plt.subplot(211)
plt.plot(vke[0].times(),vke[0].data,vae[0].times(),vae[0].data)
plt.ylabel('Vel (m/s)')
plt.legend(['Kal','Accel'])
plt.grid()
plt.xlim([0,150])
plt.subplot(212)
plt.plot(dke[0].times(),dke[0].data,dae[0].times(),dae[0].data)
plt.ylabel('Disp (m)')
plt.xlabel('Seconds')
plt.grid()
plt.xlim([0,150])
plt.suptitle('Filter corner = '+str(low_corner)+'Hz')
plt.show()

plt.figure()
plt.subplot(211)
plt.semilogx(fvke,Pvke,fvae,Pvae,fve,Pve)
plt.ylabel('Velocity PSD (dB)')
plt.legend(['Kal','Acc','GPS'])
plt.grid()
plt.subplot(212)
plt.semilogx(fdke,Pdke,fdae,Pdae,fe,Pe)
plt.ylabel('Displacement PSD (dB)')
plt.grid()
plt.xlabel('Frequency (Hz)')
plt.show()

