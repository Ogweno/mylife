from mudpy.forward import lowpass
import bbptools
from matplotlib import pyplot as plt
from numpy import genfromtxt,where

sta='st15'

stations=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/fake_nocal.stl',usecols=2,dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/run/stl/fake_nocal.stl',usecols=[0,1])
i=where(stations==sta)[0]


def ellip_lowpass(data,fcorner,fsample,order,rp=1,rs=11):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import ellip,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = ellip(order,rp,rs,array(fcorner)/(fnyquist),btype='low')
    data_filt=filtfilt(b,a,data)
    return data_filt




tlf_raw,lf_raw=bbptools.read_bbp_seismogram('/Users/dmelgar/code/BBP/bbp/bbp_data/finished/rawdata/fake_nocal/111.'+sta+'-lf.acc.000.prefilter')
tlf_bbp_proc,lf_bbp_proc=bbptools.read_bbp_seismogram('/Users/dmelgar/code/BBP/bbp/bbp_data/finished/rawdata/fake_nocal/111.'+sta+'-lf-resamp.000')
lf_filt2=lowpass(lf_raw,1.2,10,2)
lf_filt4=lowpass(lf_raw,1.2,10,4)
lf_filt8=lowpass(lf_raw,1.2,10,8)
#lf_filt2=ellip_lowpass(lf_raw,1.2,10,2,rs=0.1)
#lf_filt4=ellip_lowpass(lf_raw,1.2,10,4,rs=0.1)
#lf_filt8=ellip_lowpass(lf_raw,1.2,10,8,rs=0.1)


xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/large_eew_m5.0_frac0.5.srf')
ptime,stime=bbptools.arrivals(hypocenter,lonlat[i,0],lonlat[i,1])

#plt.figure()
#plt.plot(tlf_raw,lf_raw,'k')
#plt.plot(tlf_raw,lf_filt,'r',lw=2)
#plt.legend(['unfilt','filt'])
#plt.grid()
#
plt.figure()
plt.plot(tlf_raw,lf_filt2,'k')
plt.plot(tlf_bbp_proc,lf_bbp_proc,'b',lw=2)
plt.legend(['unfilt','bbp filt'])
plt.scatter(ptime,0,marker='|',s=80,lw=4,color='m')
plt.grid()

plt.figure()
plt.plot(tlf_raw,lf_filt2,'r')
plt.plot(tlf_raw,lf_filt4,'b')
plt.plot(tlf_raw,lf_filt8,'k')
plt.legend(['2 pole','4 pole','8 pole'])
plt.scatter(ptime,0,marker='|',s=80,lw=4,color='m')
plt.grid()

plt.show()