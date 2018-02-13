from mudpy.forward import lowpass
import bbptools
from matplotlib import pyplot as plt
from numpy import genfromtxt,where

sta='st04'
xl=[0,35]

stations=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/pwave.stl',usecols=2,dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/pwave.stl',usecols=[0,1])
i=where(stations==sta)[0]

tlf,lfe=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-lf-resamp.090')
thf,hfe=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-hf-resamp.090')
tlf,lfn=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-lf-resamp.000')
thf,hfn=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-hf-resamp.000')
tlf,lfu=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-lf-resamp.ver')
thf,hfu=bbptools.read_bbp_seismogram(u'/Users/dmelgar/code/BBP/bbp/bbp_data/tmpdata/65/65.'+sta+'-hf-resamp.ver')

xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/65/large_eew_m6.5_frac0.5.srf')
ptime,stime=bbptools.arrivals(hypocenter,lonlat[i,0],lonlat[i,1])

plt.subplot(311)
plt.plot(thf,hfn,'k')
plt.plot(tlf,lfn,'r',lw=2)
plt.legend(['HF','LF'],frameon=False)
plt.scatter(ptime,0,marker='|',s=400,lw=2)
plt.scatter(stime,0,marker='|',s=400,lw=2)
plt.ylabel('North (cm/s/s)')
plt.xlim(xl)
plt.title('M6.5, station '+sta)

plt.subplot(312)
plt.plot(thf,hfe,'k')
plt.plot(tlf,lfe,'r',lw=2)
plt.scatter(ptime,0,marker='|',s=400,lw=2)
plt.scatter(stime,0,marker='|',s=400,lw=2)
plt.ylabel('East (cm/s/s)')
plt.xlim(xl)

plt.subplot(313)
plt.plot(thf,hfu,'k')
plt.plot(tlf,lfu,'r',lw=2)
plt.scatter(ptime,0,marker='|',s=400,lw=2)
plt.scatter(stime,0,marker='|',s=400,lw=2)
plt.ylabel('Up (cm/s/s)')
plt.xlabel('Seconds after OT')
plt.xlim(xl)

plt.show()