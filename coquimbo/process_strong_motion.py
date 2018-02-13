from obspy import read
from numpy import genfromtxt,mean
from datetime import timedelta
from scipy.integrate import cumtrapz
from mudpy.forward import lowpass
from scipy.signal import filtfilt,butter
from matplotlib import pyplot as plt
from shutil import copy

#REMEMBER TO CHANGE OUTPUT PATHS FOR 0.5 VS. 1.0HZ
trimaccel=True
baseline=True
makeplots=True
move_mudpy=True

stations=genfromtxt('/Users/dmelgar/Slip_inv/Coquimbo/data/station_info/sm.gflist',usecols=0,dtype='S')
sacpath='/Users/dmelgar/Coquimbo2015/strong_motion/'
tmin=timedelta(seconds=10)
tmax=timedelta(seconds=180)
fcorner=[1./50,1./5]


#Some functions to be sued later
def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout

#Trim
if trimaccel:
    for k in range(len(stations)):
        print stations[k]
        n=read(sacpath+'proc/'+stations[k]+'.HNN.sac')
        e=read(sacpath+'proc/'+stations[k]+'.HNE.sac')
        z=read(sacpath+'proc/'+stations[k]+'.HNZ.sac')
        #Get start times
        t1=z[0].stats.starttime
        tp=z[0].stats.sac['a'] #p-wave pick on vertical
        n[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        e[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        z[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        n.write(sacpath+'trim/'+stations[k]+'.HNN.sac',format='SAC')
        e.write(sacpath+'trim/'+stations[k]+'.HNE.sac',format='SAC')
        z.write(sacpath+'trim/'+stations[k]+'.HNZ.sac',format='SAC')

#Remove pre-event mean, integrate to velocity, bandpass filter, decimate       
if baseline:
    for k in range(len(stations)):
        print stations[k]
        n=read(sacpath+'trim/'+stations[k]+'.HNN.sac')
        e=read(sacpath+'trim/'+stations[k]+'.HNE.sac')
        z=read(sacpath+'trim/'+stations[k]+'.HNZ.sac')
        #Remove zero baseline
        n[0].data=n[0].data-mean(n[0].data[0:500])
        e[0].data=e[0].data-mean(e[0].data[0:500])
        z[0].data=z[0].data-mean(z[0].data[0:500])
        #Integrate
        n[0].data=cumtrapz(n[0].data,n[0].times(),initial=0)
        e[0].data=cumtrapz(e[0].data,e[0].times(),initial=0)
        z[0].data=cumtrapz(z[0].data,z[0].times(),initial=0)
        #Lowpass or Bandpass
        if stations[k] !='CO02':
            n[0].data=lowpass(n[0].data,fcorner,1./n[0].stats.delta,2)
            e[0].data=lowpass(e[0].data,fcorner,1./e[0].stats.delta,2)
            z[0].data=lowpass(z[0].data,fcorner,1./z[0].stats.delta,2)
        else:
            n[0].data=lowpass(n[0].data,1./5,1./n[0].stats.delta,2)
            e[0].data=lowpass(e[0].data,1./5,1./e[0].stats.delta,2)
            z[0].data=lowpass(z[0].data,1./5,1./z[0].stats.delta,2)
        #Decimate to 1 Hz
        if n[0].stats.delta==0.005:
            n=stdecimate(n,5,10) ; n=stdecimate(n,5,10) ; n=stdecimate(n,4,10) ; n=stdecimate(n,2,10)
            e=stdecimate(e,5,10) ; e=stdecimate(e,5,10) ; e=stdecimate(e,4,10) ; e=stdecimate(e,2,10)
            z=stdecimate(z,5,10) ; z=stdecimate(z,5,10) ; z=stdecimate(z,4,10) ; z=stdecimate(z,2,10)
        if n[0].stats.delta==0.01:
            n=stdecimate(n,5,10) ; n=stdecimate(n,5,10) ; n=stdecimate(n,4,10)
            e=stdecimate(e,5,10) ; e=stdecimate(e,5,10) ; e=stdecimate(e,4,10)
            z=stdecimate(z,5,10) ; z=stdecimate(z,5,10) ; z=stdecimate(z,4,10)
        n.write(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNN.sac',format='SAC')
        e.write(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNE.sac',format='SAC')
        z.write(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNZ.sac',format='SAC')
        
if makeplots:
    for k in range(len(stations)):
        print stations[k]
        n=read(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNN.sac')
        e=read(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNE.sac')
        z=read(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNZ.sac')
        #n=read(sacpath+stations[k]+'.HNN.sac')
        #e=read(sacpath+stations[k]+'.HNE.sac')
        #z=read(sacpath+stations[k]+'.HNZ.sac')
        plt.figure()
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'r')
        plt.grid()
        plt.ylabel('North (m/s)')
        plt.title('Station '+stations[k]+', filter corner is 0.5Hz')
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'b')
        plt.grid()
        plt.ylabel('East (m/s)')
        plt.subplot(313)
        plt.plot(z[0].times(),z[0].data,'g')
        plt.grid()
        plt.ylabel('Up (m/s)')
        plt.savefig(sacpath+'plots/'+stations[k]+'.png')
        
if move_mudpy:
    for k in range(len(stations)):
        copy(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNN.sac','/Users/dmelgar/Slip_inv/Coquimbo_4s/data/waveforms/'+stations[k]+'.vel.n')
        copy(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNE.sac','/Users/dmelgar/Slip_inv/Coquimbo_4s/data/waveforms/'+stations[k]+'.vel.e')
        copy(sacpath+'1Hz/bandpass_0.02-0.2/'+stations[k]+'.HNZ.sac','/Users/dmelgar/Slip_inv/Coquimbo_4s/data/waveforms/'+stations[k]+'.vel.u')
    
    
    




