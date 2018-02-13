import bbptools
from numpy import genfromtxt,zeros,r_
from obspy import Stream,Trace,read
from obspy.core import UTCDateTime
import bbptools
from datetime import timedelta

make_chan_file=True
make_mseed=True

gain=1e6
chanfile='/Users/dmelgar/Cascadia_M9/cascadia_M9.chan'
data_path='/Users/dmelgar/Cascadia_M9/Cascadia_M9/bb/'
path_out='/Users/dmelgar/Cascadia_M9/planar_mseed/'
prefix='1195255'
dt_synth=0.02
noise_length=60 #in seconds
t0=UTCDateTime('2016-09-07T07:00:00')-timedelta(seconds=noise_length)

sta=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=2,dtype='S')
sta_lonlat=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=[0,1])

if make_chan_file:
    f=open(chanfile,'w')
    f.write('# net,sta,loc,chan,lat,lon,elev,samplerate,gain,units\n')
    for k in range(len(sta)):
        line1='AA\t%4s\t00\tHNE\t%.4f\t%.4f\t0.0\t50.00\t%e\tDU/M/S**2\n' % (sta[k],sta_lonlat[k,1],sta_lonlat[k,0],gain)
        line2='AA\t%4s\t00\tHNN\t%.4f\t%.4f\t0.0\t50.00\t%e\tDU/M/S**2\n' % (sta[k],sta_lonlat[k,1],sta_lonlat[k,0],gain)
        line3='AA\t%4s\t00\tHNZ\t%.4f\t%.4f\t0.0\t50.00\t%e\tDU/M/S**2\n' % (sta[k],sta_lonlat[k,1],sta_lonlat[k,0],gain)
        f.write(line1)
        f.write(line2)
        f.write(line3)
        
    f.close()
    
    
    
if make_mseed:
    noise=zeros(int(noise_length/dt_synth))
    for k in range(len(sta)):
        
        st=read(data_path+sta[k]+'.HN_.mseed')
        print 'Working on station '+sta[k]
        
        #Apply gain
        st[0].data=st[0].data*gain
        st[1].data=st[1].data*gain
        st[2].data=st[2].data*gain
        
        #Add location and network
        st[0].stats.location='00'
        st[0].stats.network='AA'   
        st[1].stats.location='00'
        st[1].stats.network='AA'   
        st[2].stats.location='00'
        st[2].stats.network='AA'        
        
        #Add noise
        st[0].data=r_[noise,st[0].data]
        st[1].data=r_[noise,st[1].data]
        st[2].data=r_[noise,st[2].data]
        
        #Change start time
        st[0].stats.starttime=t0
        st[1].stats.starttime=t0
        st[2].stats.starttime=t0
        
        #Change data type
        st[0].data = (st[0].data).astype('i4')
        st[1].data = (st[1].data).astype('i4')
        st[2].data = (st[2].data).astype('i4')
        
        #write to file
        st[0].write(path_out+sta[k]+'.'+st[0].stats.channel+'.mseed',format='MSEED',encoding=11)
        st[1].write(path_out+sta[k]+'.'+st[1].stats.channel+'.mseed',format='MSEED',encoding=11)
        st[2].write(path_out+sta[k]+'.'+st[2].stats.channel+'.mseed',format='MSEED',encoding=11)
        
