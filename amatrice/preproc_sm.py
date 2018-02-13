from obspy import Stream,Trace
from numpy import genfromtxt,r_,array
from glob import glob
from obspy.core import UTCDateTime
from obspy import read
from mudpy.forward import highpass
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt
from datetime import datetime
from string import replace

stafile='/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta'
fsta=open(stafile,'w')
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/'
fcorner=[1./20,0.4]
time_epi=UTCDateTime('2016-08-24T01:36:31.5')
tcut=25
tprevious=5
v_or_d='v'

makesac=False
signalprocess=True
makeplots_raw=False
makemudpy=True
makeplots_proc=False

def bandpass_filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'bandpass')
    data_filt=filtfilt(b,a,data)
    return data_filt
    

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


count=1
if makesac:
    files=glob(rootpath+'raw/*.ASC')
    #files=glob(rootpath+'raw/MN*.ASC')
    for k in range(len(files)):
        
        accel=genfromtxt(files[k],skip_header=64)
        accel=accel/100. #Convert to m/s/s
        st=Stream(Trace())
        
        station=files[k].split('/')[-1].split('.')[1]
        
        chan=files[k].split('/')[-1].split('.')[3]
        print 'chan before '+chan
        if chan=='HGZ':
            chan='HNZ'
        if chan=='HGE':
            chan='HNE'
        if chan=='HGN':
            chan='HNN' 
            
        if chan=='HLZ':
            chan='HNZ'
        if chan=='HLE':
            chan='HNE'
        if chan=='HLN':
            chan='HNN'  
        print '   chan after '+chan         

        st[0].data=accel
        st[0].stats.station=station
        st[0].stats.channel=chan
        
        #Read sampling interval and time of first sample
        f=open(files[k],'r')
        while True:
            line=f.readline()
            if 'DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS' in line:
                t1_string=line.split(':')[1].strip().split('.')[0]
                t1=datetime.strptime(t1_string, '%Y%m%d_%H%M%S')
            if 'STATION_LATITUDE_DEGREE' in line:
                lat=float(line.split(':')[-1])
            if 'STATION_LONGITUDE_DEGREE' in line:
                lon=float(line.split(':')[-1])
            if 'SAMPLING_INTERVAL' in line:
                dt=float(line.split(':')[1])
                break
        f.close()
        
        
        #Station file
        if chan=='HNZ':
            print count
            print station
            count+=1
            line='%s\t%.6f\t%.6f\n' % (station,lon,lat)
            fsta.write(line)

        st[0].stats.delta=dt
        st[0].stats.starttime=t1
        st.write(rootpath+'sac/'+station+'.'+chan+'.sac',format='SAC')
    fsta.close()
        
if signalprocess:
    files=glob(rootpath+'sac_w_picks/*')
    for k in range(len(files)):
        
        st=read(files[k])
        accel=st[0].data
        accel=accel-accel.mean()
        
        chan=st[0].stats.channel
        if chan=='HLZ':
            st[0].stats.channel='HNZ'
        if chan=='HLE':
            st[0].stats.channel='HNE'
        if chan=='HLN':
            st[0].stats.channel='HNN' 
        
        if st[0].stats.channel=='HNE':
            chan='e'
        if st[0].stats.channel=='HNN':
            chan='n'
        if st[0].stats.channel=='HNZ':
            chan='u'
        
        v=cumtrapz(accel,st[0].times(),initial=0)
        if v_or_d=='d':
            d=cumtrapz(v,st[0].times(),initial=0)
            fil=bandpass_filter(d,fcorner,1./st[0].stats.delta,2)
            fil=fil-fil.mean()
            st[0].data=fil
            fout=rootpath+'proc/'+st[0].stats.station+'.disp.'+chan
        else:
            fil=bandpass_filter(v,fcorner,1./st[0].stats.delta,2)
            fil=fil-fil.mean()
            st[0].data=fil
            fout=rootpath+'proc/'+st[0].stats.station+'.vel.'+chan
        
        st.write(fout,format='SAC')

if makemudpy:
    if v_or_d=='d':
        files=glob(rootpath+'proc/*disp*')
    else:
        files=glob(rootpath+'proc/*vel*')  
    for k in range(len(files)):
        st=read(files[k])
        
        #if st[0].stats.delta==0.005:
        #    st=stdecimate(st,5,2) ; st=stdecimate(st,5,2) ; st=stdecimate(st,2,2)
        #elif st[0].stats.delta==0.01:
        #    st=stdecimate(st,5,2) ; st=stdecimate(st,5,2)
        #else:
        #    print 'BURP'
        if st[0].stats.delta==0.005:
            st=stdecimate(st,5,2) ; st=stdecimate(st,4,2) ; st=stdecimate(st,2,2)
        elif st[0].stats.delta==0.01:
            st=stdecimate(st,5,2) ; st=stdecimate(st,4,2)
        else:
            print 'BURP'


        chan=st[0].stats.channel
        if chan=='HLZ':
            st[0].stats.channel='HNZ'
        if chan=='HLE':
            st[0].stats.channel='HNE'
        if chan=='HLN':
            st[0].stats.channel='HNN' 

        if st[0].stats.channel=='HNE':
            chan='e'
            reference=replace(files[k],'.e','.u')
        if st[0].stats.channel=='HNN':
            chan='n'
            reference=replace(files[k],'.n','.u')
        if st[0].stats.channel=='HNZ':
            chan='u'
            reference=files[k]
            
        st_ref=read(reference)
        ptime=st_ref[0].stats['sac']['a']
        
        #Trim them to 60s
        t1=st[0].stats.starttime+ptime-tprevious
        t2=st[0].stats.starttime+ptime+tcut
        st[0].trim(starttime=t1,endtime=t2)
        
        if v_or_d=='d':
            fout=rootpath+'mudpy/'+st[0].stats.station+'.disp.'+chan
        else:
            fout=rootpath+'mudpy/'+st[0].stats.station+'.vel.'+chan
        st.write(fout,format='SAC')        
                        
if makeplots_raw:
    stations=glob(rootpath+'sac/*HNE*')
    for k in range(len(stations)):
        sta=stations[k].split('/')[-1].split('.')[0]
        
        e=read(rootpath+'sac/'+sta+'.HNE.sac')
        n=read(rootpath+'sac/'+sta+'.HNN.sac')
        u=read(rootpath+'sac/'+sta+'.HNZ.sac')
        outname=rootpath+'plots/'+sta+'.raw.png'

        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'k')
        plt.ylabel('North m/s/s')
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'r')
        plt.ylabel('East m/s/s/')
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(313)
        plt.plot(u[0].times(),u[0].data,'m')
        plt.ylabel('Up m/s/s')
        plt.xlabel('Seconds')
        
        plt.subplots_adjust(left=0.2,right=0.97,top=0.97,bottom=0.12,wspace=0.09,hspace=0.09)
        
        plt.savefig(outname)
        plt.close()
 
    
          
if makeplots_proc:
    stations=glob(rootpath+'mudpy/*.e')
    for k in range(len(stations)):
        sta=stations[k].split('/')[-1].split('.')[0]
        
        if v_or_d=='d':
            e=read(rootpath+'mudpy/'+sta+'.disp.e')
            n=read(rootpath+'mudpy/'+sta+'.disp.n')
            u=read(rootpath+'mudpy/'+sta+'.disp.u')
            outname=rootpath+'plots/'+sta+'.proc.png'
        else:
            e=read(rootpath+'mudpy/'+sta+'.vel.e')
            n=read(rootpath+'mudpy/'+sta+'.vel.n')
            u=read(rootpath+'mudpy/'+sta+'.vel.u')
            outname=rootpath+'plots/'+sta+'.proc.vel.png'

        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'k')
        plt.ylabel('North m/s/s')
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'r')
        plt.ylabel('East m/s/s/')
        ax=plt.gca()
        ax.tick_params(labelbottom='off') 
        
        plt.subplot(313)
        plt.plot(u[0].times(),u[0].data,'m')
        plt.ylabel('Up m/s/s')
        plt.xlabel('Seconds')
        
        plt.subplots_adjust(left=0.2,right=0.97,top=0.97,bottom=0.12,wspace=0.09,hspace=0.09)
        
        plt.savefig(outname)
        plt.close()    
        
  
