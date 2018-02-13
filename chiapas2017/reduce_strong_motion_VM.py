from obspy import read
from scipy.signal import filtfilt,butter
from numpy import array,r_,zeros,genfromtxt
from glob import glob
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
    

fcorner=array([1./50,49])
    
    
stations=genfromtxt('/Users/dmelgar/Chiapas2017/strong_motion/station_info/VM_stations.txt',usecols=1,dtype='S')
tstart=genfromtxt('/Users/dmelgar/Chiapas2017/strong_motion/station_info/VM_stations.txt',usecols=5)

for k in range(len(stations)):
    f=glob('/Users/dmelgar/Chiapas2017/strong_motion/raw/_luisCM/ACEL-8Sep/*'+stations[k]+'*NE*')
    ae=read(f[0])
    f=glob('/Users/dmelgar/Chiapas2017/strong_motion/raw/_luisCM/ACEL-8Sep/*'+stations[k]+'*NN*')
    an=read(f[0])
    f=glob('/Users/dmelgar/Chiapas2017/strong_motion/raw/_luisCM/ACEL-8Sep/*'+stations[k]+'*NZ*')
    az=read(f[0])
    
    #filter
    fsample=1./an[0].stats.delta
    print fsample
    an[0].data=bandpass(an[0].data,fcorner,fsample,2)
    ae[0].data=bandpass(ae[0].data,fcorner,fsample,2)
    az[0].data=bandpass(az[0].data,fcorner,fsample,2)
    
    #de mean
    an[0].data=an[0].data-an[0].data.mean()
    ae[0].data=ae[0].data-ae[0].data.mean()
    az[0].data=az[0].data-az[0].data.mean()
    
    #To m/s/s
    an[0].data=an[0].data/100.
    ae[0].data=ae[0].data/100.
    az[0].data=az[0].data/100.
    
    #trim
    an[0].trim(starttime=an[0].stats.starttime+tstart[k],endtime=an[0].stats.starttime+tstart[k]+500)
    ae[0].trim(starttime=ae[0].stats.starttime+tstart[k],endtime=ae[0].stats.starttime+tstart[k]+500)
    az[0].trim(starttime=az[0].stats.starttime+tstart[k],endtime=az[0].stats.starttime+tstart[k]+500)
    

    
    
    an.write('/Users/dmelgar/Chiapas2017/strong_motion/sac/'+stations[k]+'.HLN.sac',format='SAC')
    ae.write('/Users/dmelgar/Chiapas2017/strong_motion/sac/'+stations[k]+'.HLE.sac',format='SAC')
    az.write('/Users/dmelgar/Chiapas2017/strong_motion/sac/'+stations[k]+'.HLZ.sac',format='SAC')
    
    plt.figure()
    plt.subplot(311)
    plt.plot(an[0].times(),an[0].data)
    plt.ylabel('North (m/s/s)')
    
    plt.subplot(312)
    plt.plot(ae[0].times(),ae[0].data)
    plt.ylabel('East (m/s/s)')
    
    plt.subplot(313)
    plt.plot(az[0].times(),az[0].data)
    plt.ylabel('Up (m/s/s)')
    plt.xlabel('Seconds')
    plt.suptitle(stations[k])
    plt.savefig('/Users/dmelgar/Chiapas2017/strong_motion/plots/'+stations[k]+'.acc.png')
    plt.close()
    
    
    