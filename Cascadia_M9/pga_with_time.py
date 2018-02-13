from obspy import read,Stream,Trace,UTCDateTime
from numpy import genfromtxt,savetxt,zeros,r_,array,c_
from run_filt import RunningMedian as med

path=u'/Users/dmelgar/Cascadia_M9/planar_mseed/'
sta=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.sta',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.sta',usecols=[1,2])
fout=u'/Users/dmelgar/code/GMT/Cascadia_M9/pga_with_time.txt'
fout2=u'/Users/dmelgar/code/GMT/Cascadia_M9/final_pga.txt'
time_epi=UTCDateTime('2016-09-07T07:00:00.000000Z')
gain=1e6

#read them all, make moving median
n=Stream()
e=Stream()
z=Stream()
for k in range(len(sta)):
    ntemp=read(path+sta[k]+'.HNN.mseed')
    ntemp.trim(starttime=time_epi)
    ntemp[0].data=(ntemp[0].data/1e6)/9.81
    etemp=read(path+sta[k]+'.HNE.mseed')
    etemp.trim(starttime=time_epi)
    etemp[0].data=(etemp[0].data/1e6)/9.81
    ztemp=read(path+sta[k]+'.HNZ.mseed')
    ztemp.trim(starttime=time_epi)
    ztemp[0].data=(ztemp[0].data/1e6)/9.81
    
    n+=ntemp[0]
    e+=etemp[0]
    z+=ztemp[0]
    
#Make time series of pga
npga=n.copy()
epga=e.copy()
zpga=z.copy()
N=n[0].stats.npts
for ksta in range(len(sta)):
    for ktime in range(N):
        if ktime==0:
            npga[ksta].data[ktime]=0
            epga[ksta].data[ktime]=0
            zpga[ksta].data[ktime]=0
        else:
            npga[ksta].data[ktime]=abs(n[ksta].data[0:ktime]).max()
            epga[ksta].data[ktime]=abs(e[ksta].data[0:ktime]).max()
            zpga[ksta].data[ktime]=abs(z[ksta].data[0:ktime]).max()
    
#Decimate time series to 2Hz
for ksta in range(len(sta)):
    npga[ksta].decimate(25,no_filter=True)
    epga[ksta].decimate(25,no_filter=True)
    zpga[ksta].decimate(25,no_filter=True)

N=npga[0].stats.npts
out=lonlat.copy()
for ktime in range(N):
    pga=zeros((len(sta),1))
    for ksta in range(len(sta)):
        pga1=npga[ksta].data[ktime]
        pga2=epga[ksta].data[ktime]
        pga3=zpga[ksta].data[ktime]
        pga[ksta]=max(pga1,pga2,pga3)

    out=c_[out,pga]
    
#Extend past the end
for k  in range(402):
    out=c_[out,pga]
    
savetxt(fout,out,fmt='%.4f')