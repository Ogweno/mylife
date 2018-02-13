from obspy import read
from numpy import mean,r_,diff
from scipy.integrate import cumtrapz

def hfilter(tr,fcorner,order):
    from scipy.signal import butter,filtfilt
    b, a = butter(order, fcorner,btype='highpass')
    y = filtfilt(b, a, tr)
    return y
    
out_path=u'/Users/dmelgar/Napa2014/acc/_baseline_v2/'
    
#gain=2.13744E+05
#n=read('/Users/dmelgar/Napa2014/acc/_raw/CE.68150..HNN.D.2014.236.102021.SAC')
#e=read('/Users/dmelgar/Napa2014/acc/_raw/CE.68150..HNE.D.2014.236.102021.SAC')
#u=read('/Users/dmelgar/Napa2014/acc/_raw/CE.68150..HNZ.D.2014.236.102021.SAC')

gain=4.26212E+05
n=read(u'/Users/dmelgar/Napa2014/acc/_raw/NC.NHC..HNE.D.2014.236.102006.SAC')
e=read(u'/Users/dmelgar/Napa2014/acc/_raw/NC.NHC..HNN.D.2014.236.102006.SAC')
u=read(u'/Users/dmelgar/Napa2014/acc/_raw/NC.NHC..HNZ.D.2014.236.102006.SAC')


n[0].data=n[0].data/gain
e[0].data=e[0].data/gain
u[0].data=u[0].data/gain


st=n.copy()
v=cumtrapz(st[0].data,st[0].times(),initial=0)
d=cumtrapz(v,st[0].times(),initial=0)
d=d-mean(d[0:300])
d=hfilter(d,1./20,2)
v=r_[0,diff(d)]/st[0].stats.delta
a=r_[0,diff(v)]/st[0].stats.delta
st[0].data=a
name=st[0].stats.station+'.'+st[0].stats.channel+'.sac'
st.write(out_path+name,format='SAC')


st=e.copy()
v=cumtrapz(st[0].data,st[0].times(),initial=0)
d=cumtrapz(v,st[0].times(),initial=0)
d=d-mean(d[0:300])
d=hfilter(d,1./20,2)
v=r_[0,diff(d)]/st[0].stats.delta
a=r_[0,diff(v)]/st[0].stats.delta
st[0].data=a
name=st[0].stats.station+'.'+st[0].stats.channel+'.sac'
st.write(out_path+name,format='SAC')


st=u.copy()
v=cumtrapz(st[0].data,st[0].times(),initial=0)
d=cumtrapz(v,st[0].times(),initial=0)
d=d-mean(d[0:300])
d=hfilter(d,1./20,2)
v=r_[0,diff(d)]/st[0].stats.delta
a=r_[0,diff(v)]/st[0].stats.delta
st[0].data=a
name=st[0].stats.station+'.'+st[0].stats.channel+'.sac'
st.write(out_path+name,format='SAC')