from obspy import Stream,Trace
from numpy import genfromtxt,r_,array
from glob import glob
from obspy.core import UTCDateTime
from obspy import read
from mudpy.forward import highpass
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt


outpath='/Users/dmelgar/Ecuador2016/strong_motion/'
t1=UTCDateTime('2016-04-16')
fcorner=[1./20,0.3]
#fcorner=[0.3]

makesac=False
signalprocess=True
integration=False
makeplots=True
makemudpy=True
plotraw=False

def filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt
    
    
def low_pass_filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'lowpass')
    data_filt=filtfilt(b,a,data)
    return data_filt

#Some functions to be sued later
def stdecimate(st,factor,order):
    from scipy.signal import butter,filtfilt
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout

if makesac:
    files=glob('/Users/dmelgar/Ecuador2016/strong_motion/raw/*')
    for k in range(len(files)):
        s=genfromtxt(files[k],skip_header=10,skip_footer=1)
        st=Stream(Trace())
        station=files[k].split('/')[-1].split('_')[0]
        print station
        chan='HN'+files[k].split('/')[-1].split('_')[2]
        accel=array([])
        for ks in range(len(s)):
            accel=r_[accel,s[ks,:]]
        st[0].data=accel
        st[0].stats.station=station
        st[0].stats.channel=chan
        st[0].stats.delta=0.01
        t1=UTCDateTime('2016-04-16')
        f=open(files[k],'r')
        while True:
            line=f.readline()
            if 'Hora' in line:
                t1.hour=int(line.split()[-3])
                t1.minute=int(line.split()[-2])
                t1.second=float(line.split()[-1])
                break
        st[0].stats.starttime=t1
        st.write(outpath+'sac/'+station+'.'+chan+'.sac',format='SAC')
        
if signalprocess:
    files=glob('/Users/dmelgar/Ecuador2016/strong_motion/sac/*')
    for k in range(len(files)):
        
        st=read(files[k])
        accel=st[0].data
        accel=accel-accel.mean()
        
        if st[0].stats.channel=='HNE':
            chan='e'
        if st[0].stats.channel=='HNN':
            chan='n'
        if st[0].stats.channel=='HNZ':
            chan='u'
        
        if integration:
            v=r_[cumtrapz(accel,st[0].times()),cumtrapz(accel,st[0].times())[-1]]
            d=r_[cumtrapz(v,st[0].times(),initial=0),cumtrapz(v,st[0].times(),initial=0)[-1]]
            fil=filter(d,fcorner,1./st[0].stats.delta,2)
            fil=fil-fil.mean()
            fil=fil/100.
            fout=outpath+'proc/'+st[0].stats.station+'.disp.'+chan
        else:
            fil=low_pass_filter(accel,fcorner,1./st[0].stats.delta,2)
            fil=fil-fil.mean()
            fil=fil/100.
            fout=outpath+'proc/'+st[0].stats.station+'.accel.'+chan
        st[0].data=fil

        
        st.write(fout,format='SAC')

if makemudpy:
    files=glob('/Users/dmelgar/Ecuador2016/strong_motion/proc/*accel*')  
    for k in range(len(files)):
        st=read(files[k])
        st=stdecimate(st,5,2) ; st=stdecimate(st,5,2) ; st=stdecimate(st,4,2)
        st[0].data=st[0].data[1:]
        if st[0].stats.channel=='HNE':
            chan='e'
        if st[0].stats.channel=='HNN':
            chan='n'
        if st[0].stats.channel=='HNZ':
            chan='u'
        #Trim them to 120s
        #t1=st[0].stats.starttime
        t1=UTCDateTime('2016-04-16T23:58:36.920')
        st[0].trim(starttime=t1,endtime=t1+120)
        
        fout=outpath+'mudpy/'+st[0].stats.station+'.accel.'+chan
        st.write(fout,format='SAC')        
                        
if makeplots:
    stations=glob('/Users/dmelgar/Ecuador2016/strong_motion/mudpy/*.e')
    for k in range(len(stations)):
        sta=stations[k].split('/')[-1].split('.')[0]
        
        plt.figure(figsize=(6,8))
        if integration:
            units='m'
            e=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.disp.e')
            n=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.disp.n')
            u=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.disp.u')
            outname=outpath+'plots/'+sta+'.disp.png'
        else:
            units='m/s/s'
            e=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.accel.e')
            n=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.accel.n')
            u=read('/Users/dmelgar/Ecuador2016/strong_motion/proc/'+sta+'.accel.u')
            outname=outpath+'plots/'+sta+'.accel.png'
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'k')
        plt.ylabel('North '+units)
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'r')
        plt.ylabel('East '+units)
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(313)
        plt.plot(u[0].times(),u[0].data,'m')
        plt.ylabel('Up '+units)
        plt.xlabel('Seconds')
        
        plt.subplots_adjust(left=0.2,right=0.97,top=0.97,bottom=0.12,wspace=0.09,hspace=0.09)
        
        plt.savefig(outname)
        plt.close()
    
if plotraw:
    files=glob('/Users/dmelgar/Ecuador2016/strong_motion/sac/*')
    for k in range(len(files)/3):
        
        plt.figure(figsize=(6,8))
        
        e=read(files[3*k])
        n=read(files[3*k+1])
        u=read(files[3*k+2])
        
        sta=files[3*k].split('/')[-1].split('.')[0]
        units='(m/s/s)'
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'k')
        plt.ylabel('North '+units)
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'r')
        plt.ylabel('East '+units)
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(313)
        plt.plot(u[0].times(),u[0].data,'m')
        plt.ylabel('Up '+units)
        plt.xlabel('Seconds')
        
        plt.subplots_adjust(left=0.2,right=0.97,top=0.97,bottom=0.12,wspace=0.09,hspace=0.09)
        
        plt.savefig(outpath+'plots/'+sta+'.raw.png')
        plt.close()     
        
  
def fixup(st,init=0,fcorner=[1./50,0.3]):
    from scipy.integrate import cumtrapz
    accel=st[0].data
    accel=accel-accel.mean()
    
    v=cumtrapz(accel,st[0].times(),initial=init)
    fil=filter(v,fcorner,1./st[0].stats.delta,2)
    fil=fil-fil.mean()
    fil=fil/100.
    
    return fil