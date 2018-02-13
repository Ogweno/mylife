'''
Spectrum compare of the 2 Kathamndu stations
'''

import nitime.algorithms as tsa
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
from numpy import log10,r_,diff,array,ones,genfromtxt
import matplotlib
from obspy.core.util.geodetics import locations2degrees,gps2DistAzimuth

matplotlib.rcParams.update({'font.size': 14})

time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=array([84.708,28.147,15]) 
kkn4_n=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXN.sac')
kkn4_e=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXE.sac')
kkn4_u=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXZ.sac')
nast_n=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXN.sac')
nast_e=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXE.sac')
nast_u=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXZ.sac')
katnp_n=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.n')
katnp_e=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.e')
katnp_u=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.vel.u')

Rf=genfromtxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/KKN4.vel.resp',usecols=0)
R=genfromtxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/KKN4.vel.resp',usecols=[1,2,3])
kkn4_Rn=R[:,0]
kkn4_Re=R[:,1]
kkn4_Ru=R[:,2]
R=genfromtxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/NAST.vel.resp',usecols=[1,2,3])
nast_Rn=R[:,0]
nast_Re=R[:,1]
nast_Ru=R[:,2]
R=genfromtxt(u'/Users/dmelgar/Nepal2015/Response_Spectra/KATNP.vel.resp',usecols=[1,2,3])
katnp_Rn=R[:,0]
katnp_Re=R[:,1]
katnp_Ru=R[:,2]


#Get distances
[d_nast,a,b]=gps2DistAzimuth(nast_n[0].stats['sac']['stla'],nast_n[0].stats['sac']['stlo'],epicenter[1],epicenter[0])
d_nast=int(d_nast/1000)
[d_kkn4,a,b]=gps2DistAzimuth(kkn4_n[0].stats['sac']['stla'],kkn4_n[0].stats['sac']['stlo'],epicenter[1],epicenter[0])
d_kkn4=int(d_kkn4/1000)
[d_katnp,a,b]=gps2DistAzimuth(katnp_n[0].stats['sac']['stla'],katnp_n[0].stats['sac']['stlo'],epicenter[1],epicenter[0])
d_katnp=int(d_katnp/1000)



kkn4_n.trim(starttime=time_epi)
kkn4_e.trim(starttime=time_epi)
kkn4_u.trim(starttime=time_epi)
nast_n.trim(starttime=time_epi)
nast_e.trim(starttime=time_epi)
nast_u.trim(starttime=time_epi)
katnp_n.trim(starttime=time_epi)
katnp_e.trim(starttime=time_epi)
katnp_u.trim(starttime=time_epi)

kkn4_n[0].data=r_[0, diff(kkn4_n[0].data)/0.2]
kkn4_e[0].data=r_[0, diff(kkn4_e[0].data)/0.2]
kkn4_u[0].data=r_[0, diff(kkn4_u[0].data)/0.2]
nast_n[0].data=r_[0, diff(nast_n[0].data)/0.2]
nast_e[0].data=r_[0, diff(nast_e[0].data)/0.2]
nast_u[0].data=r_[0, diff(nast_u[0].data)/0.2]

fn_kkn4, npsd_kkn4, nu = tsa.multi_taper_psd(kkn4_n[0].data,Fs=1./kkn4_n[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fe_kkn4, epsd_kkn4, nu = tsa.multi_taper_psd(kkn4_e[0].data,Fs=1./kkn4_e[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fu_kkn4, upsd_kkn4, nu = tsa.multi_taper_psd(kkn4_u[0].data,Fs=1./kkn4_u[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)

fn_nast, npsd_nast, nu = tsa.multi_taper_psd(nast_n[0].data,Fs=1./nast_n[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fe_nast, epsd_nast, nu = tsa.multi_taper_psd(nast_e[0].data,Fs=1./nast_e[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fu_nast, upsd_nast, nu = tsa.multi_taper_psd(nast_u[0].data,Fs=1./nast_u[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)

fn_katnp, npsd_katnp, nu = tsa.multi_taper_psd(katnp_n[0].data,Fs=1./katnp_n[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fe_katnp, epsd_katnp, nu = tsa.multi_taper_psd(katnp_e[0].data,Fs=1./katnp_e[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
fu_katnp, upsd_katnp, nu = tsa.multi_taper_psd(katnp_u[0].data,Fs=1./katnp_u[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)

#PSD in dB
#npsd_nast=10*log10(npsd_nast)
#epsd_nast=10*log10(epsd_nast)
#upsd_nast=10*log10(upsd_nast)
#npsd_kkn4=10*log10(npsd_kkn4)
#epsd_kkn4=10*log10(epsd_kkn4)
#upsd_kkn4=10*log10(upsd_kkn4)
#npsd_katnp=10*log10(npsd_katnp)
#epsd_katnp=10*log10(epsd_katnp)
#upsd_katnp=10*log10(upsd_katnp)

#Amplitude
npsd_nast=(npsd_nast)**0.5
epsd_nast=(epsd_nast)**0.5
upsd_nast=(upsd_nast)**0.5
npsd_kkn4=(npsd_kkn4)**0.5
epsd_kkn4=(epsd_kkn4)**0.5
upsd_kkn4=(upsd_kkn4)**0.5
npsd_katnp=(npsd_katnp)**0.5
epsd_katnp=(epsd_katnp)**0.5
upsd_katnp=(upsd_katnp)**0.5

#Cross correlation for lag on KATNP


fig, axarr = plt.subplots(3, 3)  

#FOr displacememnt
#ax=axarr[0,0]
#ax.plot(kkn4_n[0].times(),kkn4_n[0].data,'k')
#ax.plot(nast_n[0].times(),nast_n[0].data,'r')
#ax.legend(['KKN4','NAST'],bbox_to_anchor=(1.6, 1.35),ncol=2,frameon=False)
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticks([-0.5,0,0.5,1.0])
#ax.set_ylabel('North (m)')
#ax=axarr[1,0]
#ax.plot(kkn4_e[0].times(),kkn4_e[0].data,'k')
#ax.plot(nast_e[0].times(),nast_e[0].data,'r')
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticks([-2,0,-1,0])
#ax.set_ylabel('East (m)')
#ax=axarr[2,0]
#ax.plot(kkn4_u[0].times(),kkn4_u[0].data,'k')
#ax.plot(nast_u[0].times(),nast_u[0].data,'r')
#ax.set_ylabel('Up(m)')
#ax.yaxis.set_ticks([-0.5,0,0.5,1,1.5])
#ax.set_xlabel('Seconds after OT')

#For velocity
ax=axarr[0,0]
ax.plot(kkn4_n[0].times(),kkn4_n[0].data+2,'#606060')
ax.plot(katnp_n[0].times(),katnp_n[0].data+1,'#DAA520')
ax.plot(nast_n[0].times(),nast_n[0].data,'#DC143C')
ax.legend(['KKN4 ('+str(d_kkn4)+'km)','KATNP ('+str(d_katnp)+'km)','NAST ('+str(d_nast)+'km)'],bbox_to_anchor=(3, 1.35),ncol=3,frameon=False,fontsize=13)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticks([0,0.5,1,1.5,2.0,2.5])
ax.set_ylabel('North (m/s)')
ax=axarr[1,0]
ax.plot(kkn4_e[0].times(),kkn4_e[0].data+2,'#606060')
ax.plot(katnp_e[0].times(),katnp_e[0].data+1,'#DAA520')
ax.plot(nast_e[0].times(),nast_e[0].data,'#DC143C')
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticks([0,0.5,1,1.5,2.0,2.5])
ax.set_ylabel('East (m/s)')
ax=axarr[2,0]
ax.plot(kkn4_u[0].times(),kkn4_u[0].data+1.2,'#606060')
ax.plot(katnp_u[0].times(),katnp_u[0].data+0.6,'#DAA520')
ax.plot(nast_u[0].times(),nast_u[0].data,'#DC143C')
ax.set_ylabel('Up(m/s)')
ax.yaxis.set_ticks([-0.5,0,0.5,1.0,1.5])
ax.set_xlabel('Seconds after OT')

#Plot PSD
#ax=axarr[0,1]
#ax.semilogx(1./fn_kkn4,npsd_kkn4,'#606060')
#ax.semilogx(1./fn_katnp,npsd_katnp,'#DAA520')
#ax.semilogx(1./fn_nast,npsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])
#ax.grid(which='both')
#ax.set_ylabel('North PSD (dB)')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax=axarr[1,1]
#ax.semilogx(1./fe_kkn4,epsd_kkn4,'#606060')
#ax.semilogx(1./fe_katnp,epsd_katnp,'#DAA520')
#ax.semilogx(1./fe_nast,epsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.xaxis.set_ticklabels([])
#ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])
#ax.grid(which='both')
#ax.set_ylabel('East PSD (dB)')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax=axarr[2,1]
#ax.semilogx(1./fu_kkn4,upsd_kkn4,'#606060')
#ax.semilogx(1./fu_katnp,upsd_katnp,'#DAA520')
#ax.semilogx(1./fu_nast,upsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.set_ylabel('Up PSD (dB)')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax.grid(which='both')
#ax.set_xlabel('Period (s)')
#ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])

##Plot amplitude on linear
#ax=axarr[0,1]
#ax.plot(1./fn_kkn4,npsd_kkn4,'#606060')
#ax.plot(1./fn_katnp,npsd_katnp,'#DAA520')
#ax.plot(1./fn_nast,npsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.xaxis.set_ticklabels([])
##ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])
#ax.grid(which='both')
#ax.set_ylabel('North')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax=axarr[1,1]
#ax.plot(1./fe_kkn4,epsd_kkn4,'#606060')
#ax.plot(1./fe_katnp,epsd_katnp,'#DAA520')
#ax.plot(1./fe_nast,epsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.xaxis.set_ticklabels([])
##ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])
#ax.grid(which='both')
#ax.set_ylabel('East')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax=axarr[2,1]
#ax.plot(1./fu_kkn4,upsd_kkn4,'#606060')
#ax.plot(1./fu_katnp,upsd_katnp,'#DAA520')
#ax.plot(1./fu_nast,upsd_nast,'#DC143C')
#ax.set_xlim([1./2.5,50])
#ax.set_ylabel('Up')
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax.grid(which='both')
#ax.set_xlabel('Period (s)')
##ax.yaxis.set_ticklabels(['','-60','','-40','','-20','','0',''])

#Plot amplitude relative to KKN4
ax=axarr[0,1]
ax.plot(1./fn_katnp,ones(fn_katnp.shape),'#606060',lw=2)
ax.semilogx(1./fn_katnp,npsd_katnp/npsd_kkn4,'#DAA520')
ax.semilogx(1./fn_nast,npsd_nast/npsd_kkn4,'#DC143C')
ax.set_xlim([1./2.5,50])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','1','2','3','4','5','6'])
ax.grid(which='both')
ax.set_ylabel('North amplification')
ax=axarr[1,1]
ax.plot(1./fn_katnp,ones(fn_katnp.shape),'#606060',lw=2)
ax.semilogx(1./fe_katnp,epsd_katnp/epsd_kkn4,'#DAA520')
ax.semilogx(1./fe_nast,epsd_nast/epsd_kkn4,'#DC143C')
ax.set_xlim([1./2.5,50])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','1','2','3','4','5','6'])
ax.grid(which='both')
ax.set_ylabel('East amplification')
ax=axarr[2,1]
ax.plot(1./fn_katnp,ones(fn_katnp.shape),'#606060',lw=2)
ax.semilogx(1./fu_katnp,upsd_katnp/upsd_kkn4,'#DAA520')
ax.semilogx(1./fu_nast,upsd_nast/upsd_kkn4,'#DC143C')
ax.set_xlim([1./2.5,50])
ax.set_ylabel('Vertical amplification')
ax.grid(which='both')
ax.set_xlabel('Period (s)')
ax.yaxis.set_ticklabels(['','','1','','2','','3',''])

ax=axarr[0,2]
ax.semilogx(1./Rf,kkn4_Rn,'#606060')
ax.semilogx(1./Rf,katnp_Rn,'#DAA520')
ax.semilogx(1./Rf,nast_Rn,'#DC143C')
ax.set_xlim([1./20,20])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','','1','','2','','3'])
ax.grid(which='both')
ax.set_ylabel('North SV (m/s)')
ax=axarr[1,2]
ax.semilogx(1./Rf,kkn4_Rn,'#606060')
ax.semilogx(1./Rf,katnp_Rn,'#DAA520')
ax.semilogx(1./Rf,nast_Rn,'#DC143C')
ax.set_xlim([1./20,20])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','','1','','2','','3'])
ax.grid(which='both')
ax.set_ylabel('East SV (m/s)')
ax=axarr[2,2]
ax.semilogx(1./Rf,kkn4_Ru,'#606060')
ax.semilogx(1./Rf,katnp_Ru,'#DAA520')
ax.semilogx(1./Rf,nast_Ru,'#DC143C')
ax.set_xlim([1./20,20])
ax.yaxis.set_ticklabels(['','1','2','3','4',''])
ax.grid(which='both')
ax.set_ylabel('Vertical SV (m/s)')
ax.set_xlabel('Period (s)')

plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.85, wspace=0.25, hspace=0.05)



plt.show()