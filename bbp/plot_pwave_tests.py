from matplotlib import pyplot as plt
from obspy import read
from numpy import diff,r_
from mudpy.forward import highpass,lowpass
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from pyproj import Geod
from scipy.integrate import cumtrapz

sta='JPSB'
lonsta=-122.3467
latsta=37.1986

#sta='LUTZ'
#lonsta=-121.8652
#latsta= 37.2869

#sta='MILP'
#lonsta=-121.8340
#latsta=37.4491

path='/Users/dmelgar/FakeQuakes/M6_validation_pwave/output/waveforms/M6.000000/'

#taup
zs=8.0
g=Geod(ellps='WGS84')
azimuth,baz,dist=g.inv(-121.753508,37.332028,lonsta,latsta)
dist_in_degs=kilometer2degrees(dist/1000.)
velmod=TauPyModel('/Users/dmelgar/FakeQuakes/M6_validation_pwave/structure/bbp_norcal.npz')
Ppaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['P','p'])
p=Ppaths[0].time
Spaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['S','s'])
s=Spaths[0].time

nlf=read(path+sta+'.LYN.sac')
elf=read(path+sta+'.LYE.sac')
zlf=read(path+sta+'.LYZ.sac')

nbb=read(path+sta+'.bb.HNN.sac')
ebb=read(path+sta+'.bb.HNE.sac')
zbb=read(path+sta+'.bb.HNZ.sac')

nhf=read(path+sta+'.HNN.sac')
ehf=read(path+sta+'.HNE.sac')
zhf=read(path+sta+'.HNZ.sac')

#Diff LF to accel
dt=nlf[0].stats.delta
nlf[0].data=r_[0,diff(r_[0,diff(nlf[0].data)/dt])/dt]
elf[0].data=r_[0,diff(r_[0,diff(elf[0].data)/dt])/dt]
zlf[0].data=r_[0,diff(r_[0,diff(zlf[0].data)/dt])/dt]

#low pass filter
nlf[0].data=lowpass(nlf[0].data,1.0,1./dt,4,zerophase=True)
elf[0].data=lowpass(elf[0].data,1.0,1./dt,4,zerophase=True)
zlf[0].data=lowpass(zlf[0].data,1.0,1./dt,4,zerophase=True)

#hig pass filter
dt=nhf[0].stats.delta
nhf[0].data=highpass(nhf[0].data,1.0,1./dt,4,zerophase=True)
ehf[0].data=highpass(ehf[0].data,1.0,1./dt,4,zerophase=True)
zhf[0].data=highpass(zhf[0].data,1.0,1./dt,4,zerophase=True)


#Plot limits
elims=[-max(abs(ebb[0].data)),max(abs(ebb[0].data))]
nlims=[-max(abs(nbb[0].data)),max(abs(nbb[0].data))]
zlims=[-max(abs(zbb[0].data)),max(abs(zbb[0].data))]
xlims=[0,60]

plt.figure(figsize=(18,9))
plt.suptitle(sta)
plt.subplot(331)
plt.plot(nlf[0].times(),nlf[0].data)
plt.ylim(nlims)
plt.xlim([0,60])
plt.title('North')
plt.ylabel('a(m/s^2)')
plt.legend(['LF'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(332)
plt.plot(elf[0].times(),elf[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.title('East')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(333)
plt.plot(zlf[0].times(),zlf[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.title('Vertical')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(334)
plt.plot(nhf[0].times(),nhf[0].data)
plt.ylim(nlims)
plt.ylabel('a(m/s^2)')
plt.legend(['HF'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(335)
plt.plot(ehf[0].times(),ehf[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(336)
plt.plot(zhf[0].times(),zhf[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(337)
plt.plot(nbb[0].times(),nbb[0].data)
plt.ylim(nlims)
plt.xlim([0,60])
plt.ylabel('a(m/s^2)')
plt.xlabel('Seconds')
plt.legend(['BB'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(338)
plt.plot(ebb[0].times(),ebb[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.xlabel('Seconds')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(339)
plt.plot(zbb[0].times(),zbb[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.xlabel('Seconds')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplots_adjust(bottom=0.1,left=0.1,right=0.95,top=0.9)
plt.show()







#Dot he displacememnt thang
nlf=read(path+sta+'.LYN.sac')
elf=read(path+sta+'.LYE.sac')
zlf=read(path+sta+'.LYZ.sac')

nbb=read(path+sta+'.bb.HNN.sac')
ebb=read(path+sta+'.bb.HNE.sac')
zbb=read(path+sta+'.bb.HNZ.sac')

nhf=read(path+sta+'.HNN.sac')
ehf=read(path+sta+'.HNE.sac')
zhf=read(path+sta+'.HNZ.sac')


#Integrate everything
nhf[0].data=cumtrapz(cumtrapz(nhf[0].data,nhf[0].times(),initial=0),nhf[0].times(),initial=0)
ehf[0].data=cumtrapz(cumtrapz(ehf[0].data,ehf[0].times(),initial=0),ehf[0].times(),initial=0)
zhf[0].data=cumtrapz(cumtrapz(zhf[0].data,zhf[0].times(),initial=0),zhf[0].times(),initial=0)

nbb[0].data=cumtrapz(cumtrapz(nbb[0].data,nbb[0].times(),initial=0),nbb[0].times(),initial=0)
ebb[0].data=cumtrapz(cumtrapz(ebb[0].data,ebb[0].times(),initial=0),ebb[0].times(),initial=0)
zbb[0].data=cumtrapz(cumtrapz(zbb[0].data,zbb[0].times(),initial=0),zbb[0].times(),initial=0)

#Filters
z=False
dt=nlf[0].stats.delta
fc=1.0
nlf[0].data=lowpass(nlf[0].data,fc,1./dt,2,zerophase=z)
elf[0].data=lowpass(elf[0].data,fc,1./dt,2,zerophase=z)
zlf[0].data=lowpass(zlf[0].data,fc,1./dt,2,zerophase=z)



fc=1.
dt=nhf[0].stats.delta
nhf[0].data=highpass(nhf[0].data,fc,1./dt,2,zerophase=z)
ehf[0].data=highpass(ehf[0].data,fc,1./dt,2,zerophase=z)
zhf[0].data=highpass(zhf[0].data,fc,1./dt,2,zerophase=z)

fc=1./13
nbb[0].data=highpass(nbb[0].data,fc,1./dt,2,zerophase=z)
ebb[0].data=highpass(ebb[0].data,fc,1./dt,2,zerophase=z)
zbb[0].data=highpass(zbb[0].data,fc,1./dt,2,zerophase=z)

dt=nlf[0].stats.delta
nlf[0].data=highpass(nlf[0].data,fc,1./dt,2,zerophase=z)
elf[0].data=highpass(elf[0].data,fc,1./dt,2,zerophase=z)
zlf[0].data=highpass(zlf[0].data,fc,1./dt,2,zerophase=z)






elims=[-max(abs(ebb[0].data)),max(abs(ebb[0].data))]
nlims=[-max(abs(nbb[0].data)),max(abs(nbb[0].data))]
zlims=[-max(abs(zbb[0].data)),max(abs(zbb[0].data))]
xlims=[0,60]

plt.figure(figsize=(18,9))
plt.suptitle(sta)
plt.subplot(331)
plt.plot(nlf[0].times(),nlf[0].data)
plt.ylim(nlims)
plt.xlim([0,60])
plt.title('North')
plt.ylabel('a(m/s^2)')
plt.legend(['LF'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(332)
plt.plot(elf[0].times(),elf[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.title('East')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(333)
plt.plot(zlf[0].times(),zlf[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.title('Vertical')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(334)
plt.plot(nhf[0].times(),nhf[0].data)
plt.ylim(nlims)
plt.xlim([0,60])
plt.ylabel('a(m/s^2)')
plt.legend(['HF'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(335)
plt.plot(ehf[0].times(),ehf[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(336)
plt.plot(zhf[0].times(),zhf[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(337)
plt.plot(nbb[0].times(),nbb[0].data)
plt.ylim(nlims)
plt.xlim([0,60])
plt.ylabel('a(m/s^2)')
plt.xlabel('Seconds')
plt.legend(['BB'],frameon=False,loc=1,fontsize=16)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(338)
plt.plot(ebb[0].times(),ebb[0].data)
plt.ylim(elims)
plt.xlim([0,60])
plt.xlabel('Seconds')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplot(339)
plt.plot(zbb[0].times(),zbb[0].data)
plt.ylim(zlims)
plt.xlim([0,60])
plt.xlabel('Seconds')
plt.plot([p,p],[-10,10],'k--',lw=0.5)
plt.plot([s,s],[-10,10],'k--',lw=0.5)

plt.subplots_adjust(bottom=0.1,left=0.1,right=0.95,top=0.9)
plt.show()
