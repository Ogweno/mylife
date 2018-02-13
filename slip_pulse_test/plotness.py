from obspy import read,UTCDateTime
from string import rjust
from numpy import sqrt,array,argmax,zeros
from matplotlib import pyplot as plt
from scipy.signal import butter,filtfilt


time_epi=UTCDateTime('2017-09-08T06:11:26')
td=180

fcorner=0.5
order=2
fsample=5.

def lowpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'lowpass')
    data_filt=filtfilt(b,a,data)
        
    return data_filt


plt.figure()
N=39
tpgd=zeros(N)
for k in range(N):
    n=read(u'/Users/dmelgar/Slip_inv/Slip_pulse/output/waveforms/M8_15s.0000/PE%s.LYZ.sac' % (rjust(str(k+1),2,'0')))
    e=read(u'/Users/dmelgar/Slip_inv/Slip_pulse/output/waveforms/M8_15s.0000/PE%s.LYZ.sac' % (rjust(str(k+1),2,'0')))
    z=read(u'/Users/dmelgar/Slip_inv/Slip_pulse/output/waveforms/M8_15s.0000/PE%s.LYZ.sac' % (rjust(str(k+1),2,'0')))
    
    n[0].data=lowpass(n[0].data,fcorner,fsample,order)
    e[0].data=lowpass(e[0].data,fcorner,fsample,order)
    z[0].data=lowpass(z[0].data,fcorner,fsample,order)
    
    d=n.copy()
    d[0].data=sqrt(n[0].data**2+e[0].data**2+z[0].data**2)
    
    d[0].trim(starttime=time_epi,endtime=d[0].stats.endtime,pad=True,fill_value=0)
    
    plt.plot(d[0].times(),d[0].data)
    
    i=argmax(d[0].data)
    tpgd[k]=d[0].times()[i]
    
plt.show()