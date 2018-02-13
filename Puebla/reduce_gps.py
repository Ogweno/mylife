from obspy import read
from scipy.signal import filtfilt,butter
from numpy import array,arange
from datetime import timedelta
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

def lowpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, fcorner/fnyquist,'lowpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt
    
    
    
 
t1=10
t2=45
fcorner=array([0.3])
    
    
stations=['UTON','IGUA','TNAT','MXAS']
paths=['/Users/dmelgar/Puebla2017/GPS/trim5Hz/','/Users/dmelgar/Puebla2017/GPS/Mexico7/trim5Hz/','/Users/dmelgar/Puebla2017/GPS/trim5Hz/','/Users/dmelgar/Puebla2017/GPS/Mexico7/trim2Hz/']

for k in range(len(stations)):
    n=read(paths[k]+stations[k].lower()+'.LXN.sac')
    e=read(paths[k]+stations[k].lower()+'.LXE.sac')
    z=read(paths[k]+stations[k].lower()+'.LXZ.sac')
    
    
    
    #filter
    fsample=1./n[0].stats.delta
    print fsample
    n[0].data=lowpass(n[0].data,fcorner,fsample,2)
    e[0].data=lowpass(n[0].data,fcorner,fsample,2)
    z[0].data=lowpass(n[0].data,fcorner,fsample,2)
    
    
    tstart=z[0].stats.starttime
    tend=tstart+timedelta(seconds=t2)
    n[0].trim(starttime=tstart,endtime=tend)
    e[0].trim(starttime=tstart,endtime=tend)
    z[0].trim(starttime=tstart,endtime=tend)
    
    #resample
    t=arange(0,n[0].times()[-1],0.25)
    fn=interp1d(n[0].times(),n[0].data)
    fe=interp1d(e[0].times(),e[0].data)
    fz=interp1d(z[0].times(),z[0].data)
    ninterp=fn(t)
    einterp=fn(t)
    zinterp=fn(t)
    n[0].data=ninterp ; n[0].stats.delta=0.25
    e[0].data=ninterp ; e[0].stats.delta=0.25
    z[0].data=ninterp ; z[0].stats.delta=0.25

    
    n.write('/Users/dmelgar/Puebla2017/GPS/mudpy/'+stations[k]+'.disp.n',format='SAC')
    e.write('/Users/dmelgar/Puebla2017/GPS/mudpy/'+stations[k]+'.disp.e',format='SAC')
    z.write('/Users/dmelgar/Puebla2017/GPS/mudpy/'+stations[k]+'.disp.z',format='SAC')
    
    plt.figure()
    plt.subplot(311)
    plt.plot(n[0].times(),n[0].data)
    plt.ylabel('North (m)')
    
    plt.subplot(312)
    plt.plot(e[0].times(),e[0].data)
    plt.ylabel('East (m)')
    
    plt.subplot(313)
    plt.plot(z[0].times(),z[0].data)
    plt.ylabel('Up (m)')
    plt.xlabel('Seconds')
    plt.suptitle(stations[k])
    plt.savefig('/Users/dmelgar/Puebla2017/GPS/mudpy/_'+stations[k]+'.png')
    plt.close()
    
    
    