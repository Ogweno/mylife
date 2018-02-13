from obspy import read,UTCDateTime
from numpy import mean,where,log,log10
from mudpy.forward import lowpass as bandpass
from mtspec import mtspec
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import matplotlib as mpl
from scipy.signal import spectrogram



stations=['GO04','CO06','CO03','CO02','VA03','ROC1','VA01','MT05','VA05']
#stations=['CO06','VA03']
time_epi=UTCDateTime('2015-09-16T22:54:33')
sim_path=u'/Volumes/Illapel/FQ/illapel/output/waveforms/illapel.000002_slip7.5/'
data_path=u'/Users/dmelgar/Coquimbo2015/strong_motion/proc/'
tcut=140
fcorner=[1./10,20]
vmin=-7 ; vmax=0
cmap_spec=plt.cm.jet

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14


for k in range(len(stations)):
    
    sta=stations[k]
    print sta

    #Read fakequakes
    nfq=read(sim_path+sta+'.bb.HNN.sac')
    efq=read(sim_path+sta+'.bb.HNE.sac')
    zfq=read(sim_path+sta+'.bb.HNZ.sac')

    #cut
    nfq.trim(starttime=time_epi,endtime=time_epi+tcut)
    efq.trim(starttime=time_epi,endtime=time_epi+tcut)
    zfq.trim(starttime=time_epi,endtime=time_epi+tcut)

    #Read real data
    n=read(data_path+sta+'.HNN.sac')
    e=read(data_path+sta+'.HNE.sac')
    z=read(data_path+sta+'.HNZ.sac')
    
    #Cut real data, remove pre-event mean
    n.trim(starttime=time_epi,endtime=time_epi+tcut)
    e.trim(starttime=time_epi,endtime=time_epi+tcut)
    z.trim(starttime=time_epi,endtime=time_epi+tcut)
    n[0].data=n[0].data-mean(n[0].data[0:1000])
    e[0].data=e[0].data-mean(e[0].data[0:1000])
    z[0].data=z[0].data-mean(z[0].data[0:1000])
    
    #band pass filter both to same frequency band
    fsample=1./n[0].stats.delta
    order=4
    n[0].data=bandpass(n[0].data,fcorner,fsample,order,zerophase=True)
    e[0].data=bandpass(e[0].data,fcorner,fsample,order,zerophase=True)
    z[0].data=bandpass(z[0].data,fcorner,fsample,order,zerophase=True)
    
    fsample=1./nfq[0].stats.delta
    nfq[0].data=bandpass(nfq[0].data,fcorner,fsample,order,zerophase=True)
    efq[0].data=bandpass(efq[0].data,fcorner,fsample,order,zerophase=True)
    zfq[0].data=bandpass(zfq[0].data,fcorner,fsample,order,zerophase=True)
    
    
    #Get psds
    Tw=3.5
    Tw_spec=3.5
    Ntapers=5
    dt=n[0].stats.delta
    N, fn, jackknife, _, _ = mtspec(data=n[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=n[0].stats.npts, statistics=True)
    E, fe, jackknife, _, _ = mtspec(data=e[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=e[0].stats.npts, statistics=True)
    Z, fz, jackknife, _, _ = mtspec(data=z[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=z[0].stats.npts, statistics=True)
    
    dt=nfq[0].stats.delta
    Nfq, ffq, jackknife, _, _ = mtspec(data=nfq[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=nfq[0].stats.npts, statistics=True)
    Efq, ffq, jackknife, _, _ = mtspec(data=efq[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=efq[0].stats.npts, statistics=True)
    Zfq, ffq, jackknife, _, _ = mtspec(data=zfq[0].data, delta=dt, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=zfq[0].stats.npts, statistics=True)
    
    #Make amplitude
    N=N**0.5 ; E=E**0.5 ; Z=Z**0.5
    Nfq=Nfq**0.5 ; Efq=Efq**0.5 ; Zfq=Zfq**0.5
    
    
    #Get spectrograms
    nperseg=256
    noverlap=128
    fsample=1./n[0].stats.delta
    ffn, ttn, Nspec = spectrogram(n[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')
    ffe, tte, Espec = spectrogram(e[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')
    ffz, ttz, Zspec = spectrogram(z[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')

    fsample=1./nfq[0].stats.delta
    ff_fq, tt_fq, Nspec_fq = spectrogram(nfq[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')
    ff_fq, tt_fq, Espec_fq = spectrogram(efq[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')
    ff_fq, tt_fq, Zspec_fq = spectrogram(zfq[0].data, fsample,nperseg=nperseg,noverlap=noverlap,scaling='spectrum')
    
    #Get rid of zero frequency
    ffn[0]=0.1
    ffe[0]=0.1
    ffz[0]=0.1
    #
    ff_fq[0]=0.1#ff_fq[1:]
    #Nspec_fq=Nspec_fq[1:,:]
    #Espec_fq=Espec_fq[1:,:]
    #Zspec_fq=Zspec_fq[1:,:]
    
    #resample real spectra to fakequakes frequencies
    interp=interp1d(fn,N)
    N=interp(ffq)
    interp=interp1d(fe,E)
    E=interp(ffq)
    interp=interp1d(fz,Z)
    Z=interp(ffq)
    
    #Keep only frequencies of itnerst
    i=where((ffq>0.1) & (ffq<20))[0]
    ffq=ffq[i]
    N=N[i]; E=E[i] ; Z=Z[i]
    Nfq=Nfq[i] ; Efq=Efq[i] ; Zfq=Zfq[i]
    
    
    
    fig, axarr = plt.subplots(6, 3,figsize=(16,16)) 
    lw=0.4
    lwf=1.2
    flims=[0.1,20]
    
    #North
    ymax=max(abs(n[0].data).max(),abs(nfq[0].data).max())
    ax=axarr[0,0]
    ax.plot(n[0].times(),n[0].data,c='#DC143C',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.set_title('North',fontsize=14)
    ax.annotate('Obs',xy=(10,0.6*ymax),fontsize=14)
    ax.set_ylabel(r'a($m/s^2$)',fontsize=14)
    ax.set_xlabel('time (s)',fontsize=14)
    
    ax=axarr[1,0]
    ax.plot(nfq[0].times(),nfq[0].data,c='#6495ED',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.annotate('FQ',xy=(10,0.6*ymax),fontsize=14)
    ax.set_xlabel('time (s)',fontsize=14)
    ax.set_ylabel(r'a($m/s^2$)',fontsize=14)
    
    ax=axarr[4,0]
    ax.loglog(ffq,N,c='#DC143C',lw=lwf)
    ax.loglog(ffq,Nfq,c='#6495ED',lw=lwf)
    ax.set_xlim(flims)
    ax.set_ylabel('Four. Ampl.',fontsize=14)
    ax.legend(['Obs','FQ'],frameon=False,loc=3,fontsize=14,ncol=2)
    ax.set_xlabel('f (Hz)',fontsize=14)
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[5,0]
    ax.semilogx(ffq,log(N/Nfq),c='#708090',lw=1)
    ax.set_xlim(flims)
    ax.set_ylabel('Ln(obs/FQ)',fontsize=14)
    ax.set_xlabel('f (Hz)',fontsize=14)
    ax.set_ylim([-3,3])
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[2,0]
    ax.set_yscale('log')
    ax.annotate('Obs',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(ttn,ffn,log10(Nspec),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylabel('f (Hz)',fontsize=14)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    ax=axarr[3,0]
    ax.set_yscale('log')
    ax.annotate('FQ',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(tt_fq,ff_fq,log10(Nspec_fq),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylabel('f (Hz)',fontsize=14)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    
    
    
    
    
    #East
    ymax=max(abs(e[0].data).max(),abs(efq[0].data).max())
    ax=axarr[0,1]
    ax.plot(e[0].times(),e[0].data,c='#DC143C',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.set_title('East',fontsize=14)
    ax.annotate('Obs',xy=(10,0.6*ymax),fontsize=14)
    ax.set_xlabel('time (s)',fontsize=14)
    
    ax=axarr[1,1]
    ax.plot(efq[0].times(),efq[0].data,c='#6495ED',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.annotate('FQ',xy=(10,0.6*ymax),fontsize=14)
    ax.set_xlabel('time (s)',fontsize=14)
    
    ax=axarr[4,1]
    ax.loglog(ffq,E,c='#DC143C',lw=lwf)
    ax.loglog(ffq,Efq,c='#6495ED',lw=lwf)
    ax.set_xlim(flims)
    ax.set_xlabel('f (Hz)',fontsize=14)
    ax.legend(['Obs','FQ'],frameon=False,loc=3,fontsize=14,ncol=2)
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[5,1]
    ax.semilogx(ffq,log(E/Efq),c='#708090',lw=1)
    ax.set_xlim(flims)
    ax.set_xlabel('f (Hz)',fontsize=14)
    ax.set_ylim([-3,3])
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[2,1]
    ax.set_yscale('log')
    ax.annotate('Obs',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(tte,ffe,log10(Espec),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    ax=axarr[3,1]
    ax.set_yscale('log')
    ax.annotate('FQ',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(tt_fq,ff_fq,log10(Espec_fq),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    
    #Up
    ymax=max(abs(z[0].data).max(),abs(zfq[0].data).max())
    ax=axarr[0,2]
    ax.plot(z[0].times(),z[0].data,c='#DC143C',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.set_title('Up',fontsize=14)
    ax.annotate('Obs',xy=(10,0.6*ymax),fontsize=14)
    ax.set_xlabel('time (s)',fontsize=14)
    
    ax=axarr[1,2]
    ax.plot(zfq[0].times(),zfq[0].data,c='#6495ED',lw=lw)
    ax.set_ylim([-ymax,ymax])
    ax.set_xlim([0,tcut])
    ax.annotate('FQ',xy=(10,0.6*ymax),fontsize=14)
    ax.set_xlabel('Seconds after OT',fontsize=14)
    
    ax=axarr[4,2]
    ax.loglog(ffq,Z,c='#DC143C',lw=lwf)
    ax.loglog(ffq,Zfq,c='#6495ED',lw=lwf)
    ax.set_xlim(flims)
    ax.set_xlabel('f (Hz)',fontsize=14)
    ax.legend(['Obs','FQ'],frameon=False,loc=3,fontsize=14,ncol=2)
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[5,2]
    ax.semilogx(ffq,log(Z/Zfq),c='#708090',lw=1)
    ax.set_xlim(flims)
    ax.set_xlabel('f (Hz)',fontsize=14)
    ax.set_ylim([-3,3])
    yl=ax.get_ylim()
    ax.plot([1,1],yl,'--',lw=2.5,c='#C0C0C0')
    
    ax=axarr[2,2]
    ax.set_yscale('log')
    ax.annotate('Obs',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(ttz,ffz,log10(Zspec),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    ax=axarr[3,2]
    ax.set_yscale('log')
    ax.annotate('FQ',xy=(4,6.58),fontsize=14,color='w')
    ax.pcolormesh(tt_fq,ff_fq,log10(Zspec_fq),cmap=cmap_spec,vmin=vmin,vmax=vmax)
    ax.set_ylim(flims)
    ax.set_xlim([0,tcut])
    ax.set_xlabel('time (s)',fontsize=14)
    ax.plot([0,tcut],[1,1],'--',lw=2.5,c='#A0A0A0')
    
    plt.suptitle('Station '+ sta,fontsize=16)
    
    fig.subplots_adjust(left=0.07,right=0.98,bottom=0.05,top=0.95,hspace=0.4)
    fig.savefig(sim_path+'_'+sta+'.png')
    plt.close()
