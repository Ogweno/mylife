from obspy import read
from scipy.signal import filtfilt,butter
from numpy import array,r_,zeros
from datetime import timedelta
from matplotlib import pyplot as plt
from scipy.integrate import cumtrapz

def bandpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, fcorner/fnyquist,'bandpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt
    
    
def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout
    
 
t1=10
t2=45
fcorner=array([1./20,0.4])
    
    
stations=['FTIG','HLIG','TOIG','TPIG','PPIG','YAIG','PLIG','ARIG','MEIG','TLIG']

for k in range(len(stations)):
    an=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/_withPicks/'+stations[k]+'.HLN.sac')
    ae=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/_withPicks/'+stations[k]+'.HLE.sac')
    az=read(u'/Users/dmelgar/Puebla2017/strong_motion/sac/_withPicks/'+stations[k]+'.HLZ.sac')
    p=az[0].stats['sac']['a']
    
    #Integrate
    n=an.copy()
    e=an.copy()
    z=an.copy()
    
    n[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    e[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    z[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    n[0].data=cumtrapz(n[0].data,n[0].times(),initial=0)
    e[0].data=cumtrapz(e[0].data,e[0].times(),initial=0)
    z[0].data=cumtrapz(z[0].data,z[0].times(),initial=0)
    
    #filter
    fsample=1./n[0].stats.delta
    print fsample
    n[0].data=bandpass(n[0].data,fcorner,fsample,2)
    e[0].data=bandpass(e[0].data,fcorner,fsample,2)
    z[0].data=bandpass(z[0].data,fcorner,fsample,2)
    
    
    tstart=z[0].stats.starttime+timedelta(seconds=int(p-t1))
    tend=tstart+timedelta(seconds=t2)
    n[0].trim(starttime=tstart,endtime=tend)
    e[0].trim(starttime=tstart,endtime=tend)
    z[0].trim(starttime=tstart,endtime=tend)
    
    #decimate
    n1=stdecimate(n,5,2) ; n2=stdecimate(n1,5,2)
    e1=stdecimate(e,5,2) ; e2=stdecimate(e1,5,2)
    z1=stdecimate(z,5,2) ; z2=stdecimate(z1,5,2)

    #delay
    if stations[k]=='TPIG':
        nsamples=int(3/0.25)
        print 'TPIG'
        n2[0].data=r_[zeros(nsamples),n2[0].data[0:-nsamples]]
        e2[0].data=r_[zeros(nsamples),e2[0].data[0:-nsamples]] 
        z2[0].data=r_[zeros(nsamples),z2[0].data[0:-nsamples]]
    
    
    n2.write('/Users/dmelgar/Puebla2017/strong_motion/mudpy/'+stations[k]+'.disp.n',format='SAC')
    e2.write('/Users/dmelgar/Puebla2017/strong_motion/mudpy/'+stations[k]+'.disp.e',format='SAC')
    z2.write('/Users/dmelgar/Puebla2017/strong_motion/mudpy/'+stations[k]+'.disp.u',format='SAC')
    
    plt.figure()
    plt.subplot(311)
    plt.plot(n2[0].times(),n2[0].data)
    plt.ylabel('North (m)')
    
    plt.subplot(312)
    plt.plot(e2[0].times(),e2[0].data)
    plt.ylabel('East (m)')
    
    plt.subplot(313)
    plt.plot(z2[0].times(),z2[0].data)
    plt.ylabel('Up (m)')
    plt.xlabel('Seconds')
    plt.suptitle(stations[k])
    plt.savefig('/Users/dmelgar/Puebla2017/strong_motion/mudpy/_'+stations[k]+'disp.png')
    plt.close()
    
    
    