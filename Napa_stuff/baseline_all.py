from glob import glob
from obspy import read
from numpy import mean,r_,diff
from scipy.integrate import cumtrapz

def hfilter(tr,fcorner,order):
    from scipy.signal import butter,filtfilt
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y
    
def lfilter(tr,fcorner,order):
    from scipy.signal import butter,filtfilt
    b, a = butter(order, fcorner,btype='lowpass')
    y = filtfilt(b, a, tr)
    return y
    
def bpfilter(tr,fcorner,order):
    from scipy.signal import butter,filtfilt
    b, a = butter(order, fcorner,btype='bandpass')
    y = filtfilt(b, a, tr)
    return y
    
    
files=glob(u'/Users/dmelgar/Napa2014/acc/_trim/*.sac')
out_path=u'/Users/dmelgar/Napa2014/acc/_baseline_v1/'

# version 1 is fdemean and filter
# version 2 is integrate, de mean, filter then diff


for k in range(len(files)):
    print files[k]
    st=read(files[k])
    #v=cumtrapz(st[0].data,st[0].times(),initial=0)
    #d=cumtrapz(v,st[0].times(),initial=0)
    #d=d-mean(d[0:300])
    #d=hfilter(d,1./20,2)
    #v=r_[0,diff(d)]/st[0].stats.delta
    #a=r_[0,diff(v)]/st[0].stats.delta
    a=st[0].data-mean(st[0].data[0:500])
    a=hfilter(a,1./20,2)
    st[0].data=a
    name=st[0].stats.station+'.'+st[0].stats.channel+'.sac'
    st.write(out_path+name,format='SAC')