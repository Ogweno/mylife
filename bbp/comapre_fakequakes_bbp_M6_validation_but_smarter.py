import bbptools
from mudpy.forward import highpass,lowpass
from obspy import read
from numpy import r_,diff,genfromtxt
from matplotlib import pyplot as plt

sta='MONB'
xl=[0,30]
run='7632866'
M='M5'
lag=0.5

#read in BBP waveforms
tlf,lfe=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-lf.acc.090')
thf,hfe=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-hf.acc.090')
tlf,lfn=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-lf.acc.000')
thf,hfn=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-hf.acc.000')
tlf,lfz=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-lf.acc.ver')
thf,hfz=bbptools.read_bbp_seismogram(u'/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/tmpdata/'+run+'/'+run+'.'+sta+'-hf.acc.ver')

#Read full waveforms
s=genfromtxt('/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/outdata/'+run+'/'+run+'.'+sta+'.acc.bbp')
tbb=s[:,0] ; nbb=s[:,1] ; ebb=s[:,2] ; zbb=s[:,3]

#filter
fsample=1/(tlf[1]-tlf[0])
lfe=lowpass(lfe,1.0,fsample,2,zerophase=True)
lfn=lowpass(lfn,1.0,fsample,2,zerophase=True)
lfz=lowpass(lfz,1.0,fsample,2,zerophase=True)

fsample=1/(thf[1]-thf[0])
hfe=highpass(hfe,1.0,fsample,2,zerophase=True)
hfn=highpass(hfn,1.0,fsample,2,zerophase=True)
hfz=highpass(hfz,1.0,fsample,2,zerophase=True)

#rescale
lfe=lfe/100
lfn=lfn/100
lfz=lfz/100
hfe=hfe/100
hfn=hfn/100
hfz=hfz/100
nbb=nbb/100
ebb=ebb/100
zbb=zbb/100


#read fakequakes

fq_lfe=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.LYE.sac')
fq_lfn=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.LYN.sac')
fq_lfz=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.LYZ.sac')
fq_hfe=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.HNE.sac')
fq_hfn=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.HNN.sac')
fq_hfz=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.HNZ.sac')
fq_bbe=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.bb.HNE.sac')
fq_bbn=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.bb.HNN.sac')
fq_bbz=read(u'/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/M6.000000/'+sta+'.bb.HNZ.sac')

# LF to acceleration
dt=fq_lfn[0].stats.delta

fq_lfn[0].data=r_[0,diff(fq_lfn[0].data)/dt]
fq_lfe[0].data=r_[0,diff(fq_lfe[0].data)/dt]
fq_lfz[0].data=r_[0,diff(fq_lfz[0].data)/dt]

fq_lfn[0].data=r_[0,diff(fq_lfn[0].data)/dt]
fq_lfe[0].data=r_[0,diff(fq_lfe[0].data)/dt]
fq_lfz[0].data=r_[0,diff(fq_lfz[0].data)/dt]

#lowpass
fsample=1/dt
fq_lfe[0].data=lowpass(fq_lfe[0].data,1.0,fsample,4,zerophase=True)
fq_lfn[0].data=lowpass(fq_lfn[0].data,1.0,fsample,4,zerophase=True)
fq_lfz[0].data=lowpass(fq_lfz[0].data,1.0,fsample,4,zerophase=True)

#highpass
fsample=1/fq_hfe[0].stats.delta
fq_hfe[0].data=highpass(fq_hfe[0].data,1.0,fsample,4,zerophase=True)
fq_hfn[0].data=highpass(fq_hfn[0].data,1.0,fsample,4,zerophase=True)
fq_hfz[0].data=highpass(fq_hfz[0].data,1.0,fsample,4,zerophase=True)




# Make plots

#LF plot
plt.figure(figsize=(9,7))

plt.subplot(311)
plt.plot(tlf,lfe,c='k')
plt.plot(fq_lfe[0].times(),fq_lfe[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('East (m/s/s)')
plt.legend(['SCEC BBP','FakeQuakes'])
plt.xlim(xl)

plt.subplot(312)
plt.plot(tlf,lfn,c='k')
plt.plot(fq_lfn[0].times(),fq_lfn[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('North (m/s/s)')
plt.xlim(xl)

plt.subplot(313)
plt.plot(tlf,lfz,c='k')
plt.plot(fq_lfz[0].times(),fq_lfz[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('Up (m/s/s)')
plt.xlim(xl)





#HF plot
plt.figure(figsize=(9,7))

plt.subplot(311)
plt.plot(thf,hfe,c='k')
plt.plot(fq_hfe[0].times()-lag,fq_hfe[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('East (m/s/s)')
plt.legend(['SCEC BBP','FakeQuakes'])
plt.xlim(xl)

plt.subplot(312)
plt.plot(thf,hfn,c='k')
plt.plot(fq_hfn[0].times()-lag,fq_hfn[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('North (m/s/s)')
plt.xlim(xl)

plt.subplot(313)
plt.plot(thf,hfz,c='k')
plt.plot(fq_hfz[0].times()-lag,fq_hfz[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('Up (m/s/s)')
plt.xlim(xl)



#BB plot
plt.figure(figsize=(9,7))

plt.subplot(311)
plt.plot(tbb,ebb,c='k')
plt.plot(fq_bbe[0].times(),fq_bbe[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('East (m/s/s)')
plt.legend(['SCEC BBP','FakeQuakes'])
plt.xlim(xl)

plt.subplot(312)
plt.plot(tbb,nbb,c='k')
plt.plot(fq_bbn[0].times(),fq_bbn[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('North (m/s/s)')
plt.xlim(xl)

plt.subplot(313)
plt.plot(tbb,zbb,c='k')
plt.plot(fq_bbz[0].times(),fq_bbz[0].data,c='r')
plt.xlabel('Seconds')
plt.ylabel('Up (m/s/s)')
plt.xlim(xl)


plt.show()