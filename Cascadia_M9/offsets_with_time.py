from obspy import read,Stream,Trace
from numpy import genfromtxt,savetxt,zeros,r_,array,c_
from run_filt import RunningMedian as med

path=u'/Users/dmelgar/Slip_inv/Cascadia_M9/output/waveforms/planar_lowpass.000006/'
sta=genfromtxt(u'/Users/dmelgar/Slip_inv/Cascadia_M9/data/station_info/cascadia_small.gflist',usecols=0,dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/Slip_inv/Cascadia_M9/data/station_info/cascadia_small.gflist',usecols=[1,2])
fout=u'/Users/dmelgar/code/GMT/Cascadia_M9/offsets_with_time.txt'

Nmed=10

#read them all, make moving median
nraw=Stream()
eraw=Stream()
zraw=Stream()
pad=zeros(Nmed-1)
for k in range(len(sta)):
    nraw+=read(path+sta[k]+'.LYN.sac')
    eraw+=read(path+sta[k]+'.LYE.sac')
    zraw+=read(path+sta[k]+'.LYZ.sac')
    
#Now get moving median
n=nraw.copy()
e=eraw.copy()
z=zraw.copy()
for k in range(len(sta)):
    n[k].data=r_[pad,array(med(nraw[k],Nmed))]
    e[k].data=r_[pad,array(med(eraw[k],Nmed))]
    z[k].data=r_[pad,array(med(zraw[k],Nmed))]
    

N=n[0].stats.npts
out=lonlat.copy()
for ktime in range(N):
    north=zeros((len(sta),1))
    east=zeros((len(sta),1))
    up=zeros((len(sta),1))
    for ksta in range(len(sta)):
        north[ksta]=n[ksta].data[ktime]
        east[ksta]=e[ksta].data[ktime]
        up[ksta]=z[ksta].data[ktime]
    
    out=c_[out,east,north,up,east,north,up]  #do it twice because gmt time slices are very 0.5s
    
savetxt(fout,out,fmt='%.4f')