from obspy import read,Stream,Trace,UTCDateTime
from numpy import genfromtxt,mean
from datetime import timedelta
from scipy.integrate import cumtrapz
from mudpy.forward import lowpass
from scipy.signal import filtfilt,butter
from matplotlib import pyplot as plt
from shutil import copy
from glob import glob

make_station_file=False
makesac=False
trimaccel=True
baseline=True
makeplots=True
move_mudpy=True

station_file='/Users/dmelgar/Melinka2016/station_info/sm.sta'
path='/Users/dmelgar/Melinka2016/strong_motion/'
tmin=timedelta(seconds=10)
tmax=timedelta(seconds=70)
fcorner=[1./50,0.3]


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

if make_station_file:
    file_list=glob(path+'raw/*HNE*')
    fout=open(station_file,'w')
    fout.write('# station,lon,lat\n')
    for k in range(len(file_list)):
        f=open(file_list[k],'r')
        while True:
            line=f.readline()
            if 'Estacion' in line:
                sta=line.split()[2]
            if 'Latitud' in line:
                lon=float(line.split()[4])
                lat=float(line.split()[2])
                fout.write('%s\t%.6f\t%.6f\n' % (sta,lon,lat))
                break
    fout.close()
    f.close()
        
if makesac:
    file_list=glob(path+'raw/*.txt')
    for k in range(len(file_list)):
        print k
        station=file_list[k].split('-')[2]
        chan=file_list[k].split('-')[-1].split('.')[0]
        f=open(file_list[k],'r')
        while True:
            line=f.readline()
            if 'Tiempo' in line:
                start_time=line.split()[-1]
            if  'Tasa' in line:
                sample_rate=float(line.split()[-2])
                break
                
        data=genfromtxt(file_list[k])
        
        st=Stream(Trace())
        st[0].data=data
        st[0].stats.station=station
        st[0].stats.channel=chan
        st[0].stats.delta=1./sample_rate
        st[0].stats.starttime=UTCDateTime(start_time)
        filename=station+'.'+chan+'.sac'
        st.write(path+'sac/'+filename,format='SAC')
    


#Trim
if trimaccel:
    stations=genfromtxt(u'/Users/dmelgar/Melinka2016/station_info/sm.sta',usecols=0,dtype='S')
    for k in range(len(stations)):
        print stations[k]
        n=read(path+'sac/'+stations[k]+'.HNN.sac')
        e=read(path+'sac/'+stations[k]+'.HNE.sac')
        z=read(path+'sac/'+stations[k]+'.HNZ.sac')
        #Get start times
        t1=z[0].stats.starttime
        tp=z[0].stats.sac['a'] #p-wave pick on vertical
        n[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        e[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        z[0].trim(starttime=t1+tp-tmin,endtime=t1+tp+tmax)
        n.write(path+'trim/'+stations[k]+'.HNN.sac',format='SAC')
        e.write(path+'trim/'+stations[k]+'.HNE.sac',format='SAC')
        z.write(path+'trim/'+stations[k]+'.HNZ.sac',format='SAC')

#Remove pre-event mean, integrate to velocity, bandpass filter, decimate       
if baseline:
    stations=genfromtxt(u'/Users/dmelgar/Melinka2016/station_info/sm.sta',usecols=0,dtype='S')
    for k in range(len(stations)):
        print stations[k]
        n=read(path+'trim/'+stations[k]+'.HNN.sac')
        e=read(path+'trim/'+stations[k]+'.HNE.sac')
        z=read(path+'trim/'+stations[k]+'.HNZ.sac')
        #Remove zero baseline
        n[0].data=n[0].data-mean(n[0].data[0:500])
        e[0].data=e[0].data-mean(e[0].data[0:500])
        z[0].data=z[0].data-mean(z[0].data[0:500])
        #Integrate
        n[0].data=cumtrapz(n[0].data,n[0].times(),initial=0)
        e[0].data=cumtrapz(e[0].data,e[0].times(),initial=0)
        z[0].data=cumtrapz(z[0].data,z[0].times(),initial=0)
        #Lowpass or Bandpass
        n[0].data=lowpass(n[0].data,fcorner,1./n[0].stats.delta,2)
        e[0].data=lowpass(e[0].data,fcorner,1./e[0].stats.delta,2)
        z[0].data=lowpass(z[0].data,fcorner,1./z[0].stats.delta,2)
        #Decimate to 4 Hz
        if n[0].stats.delta==0.005:
            n=stdecimate(n,5,10) ; n=stdecimate(n,5,10) ; n=stdecimate(n,2,10)
            e=stdecimate(e,5,10) ; e=stdecimate(e,5,10) ; e=stdecimate(e,2,10) 
            z=stdecimate(z,5,10) ; z=stdecimate(z,5,10) ; z=stdecimate(z,2,10) 
        if n[0].stats.delta==0.01:
            n=stdecimate(n,5,10) ; n=stdecimate(n,5,10)
            e=stdecimate(e,5,10) ; e=stdecimate(e,5,10)
            z=stdecimate(z,5,10) ; z=stdecimate(z,5,10)
        n.write(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNN.sac',format='SAC')
        e.write(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNE.sac',format='SAC')
        z.write(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNZ.sac',format='SAC')
        
if makeplots:
    stations=genfromtxt(u'/Users/dmelgar/Melinka2016/station_info/sm.sta',usecols=0,dtype='S')
    for k in range(len(stations)):
        print stations[k]
        n=read(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNN.sac')
        e=read(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNE.sac')
        z=read(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNZ.sac')
        #n=read(path+'sac/'+stations[k]+'.HNN.sac')
        #e=read(path+'sac/'+stations[k]+'.HNE.sac')
        #z=read(path+'sac/'+stations[k]+'.HNZ.sac')
        plt.figure()
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'r')
        plt.grid()
        plt.ylabel('North (m/s/s)')
        plt.title('Station '+stations[k]+', filter corners are 0.02-0.3Hz')
        #plt.title('Station '+stations[k]+', no filter')
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'b')
        plt.grid()
        plt.ylabel('East (m/s/s)')
        plt.subplot(313)
        plt.plot(z[0].times(),z[0].data,'g')
        plt.grid()
        plt.ylabel('Up (m/s/s)')
        plt.savefig(path+'plots/'+stations[k]+'.baseline.png')
        plt.close()
        
if move_mudpy:
    for k in range(len(stations)):
        copy(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNN.sac','/Users/dmelgar/Slip_inv/Melinka/data/waveforms/'+stations[k]+'.vel.n')
        copy(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNE.sac','/Users/dmelgar/Slip_inv/Melinka/data/waveforms/'+stations[k]+'.vel.e')
        copy(path+'4Hz/bandpass_0.02-0.3/'+stations[k]+'.HNZ.sac','/Users/dmelgar/Slip_inv/Melinka/data/waveforms/'+stations[k]+'.vel.u')
    
    
    




