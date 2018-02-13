from bbptools import paste_seismograms
from numpy import genfromtxt
from obspy import Stream,Trace
from obspy.core import UTCDateTime

#Paths
lf_path=u'/Users/dmelgar/Cascadia_M9/Cascadia_M9/lf/'
hf_path=u'/Users/dmelgar/Cascadia_M9/Cascadia_M9/hf/'
out_path=u'/Users/dmelgar/Cascadia_M9/Cascadia_M9/bb/'

#Start time
tstart=UTCDateTime('2016-09-07T14:42:00')

#load station names
sta=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=2,dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=[0,1])

for k in range(len(sta)):
    
    print sta[k]
    
    # VELOCITY
    

    
    
    # Accel
    
    #east
    lf_file=lf_path+sta[k]+'.090'
    hf_file=hf_path+'4367304.'+sta[k]+'-hf.090'
    t,seis=paste_seismograms(lf_file,hf_file,filter_corner=1.0,filter_order=2,two_pass=False,seismogram_type='a',dt_synth=0.02)
    seis=seis/100.
    out=Stream(Trace())
    out[0].data=seis
    out[0].stats.delta=t[1]-t[0]
    out[0].stats.starttime=tstart
    out[0].stats.station=sta[k]
    out[0].stats.channel='HNE'
    
    #north
    lf_file=lf_path+sta[k]+'.000'
    hf_file=hf_path+'4367304.'+sta[k]+'-hf.000'
    t,seis=paste_seismograms(lf_file,hf_file,filter_corner=1.0,filter_order=2,two_pass=False,seismogram_type='a',dt_synth=0.02)
    seis=seis/100.
    out+=Trace()
    out[1].data=seis
    out[1].stats.delta=t[1]-t[0]
    out[1].stats.starttime=tstart
    out[1].stats.station=sta[k]
    out[1].stats.channel='HNN'
    
    #up
    lf_file=lf_path+sta[k]+'.ver'
    hf_file=hf_path+'4367304.'+sta[k]+'-hf.ver'
    t,seis=paste_seismograms(lf_file,hf_file,filter_corner=1.0,filter_order=2,two_pass=False,seismogram_type='a',dt_synth=0.02)
    seis=seis/100.
    out+=Trace()
    out[2].data=seis
    out[2].stats.delta=t[1]-t[0]
    out[2].stats.starttime=tstart
    out[2].stats.station=sta[k]
    out[2].stats.channel='HNZ'
    out.write(out_path+sta[k]+'.HN_.mseed',format='MSEED')
    
    
 
    