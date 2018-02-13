from obspy import read,Stream,Trace,UTCDateTime
from numpy import genfromtxt,savetxt,zeros,r_,array,c_
from run_filt import RunningMedian as med

path=u'/Users/dmelgar/Cascadia_M9/planar_mseed/'
sta=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.sta',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Cascadia_M9/pnw.planar.0006.sta',usecols=[1,2])
fout=u'/Users/dmelgar/Cascadia_M9/pga_final.txt'
time_epi=UTCDateTime('2016-09-07T07:00:00.000000Z')
gain=1.e6

f=open(fout,'w')
f.write('# sta, lon, lat, PGA (g)\n')
#read them all
for k in range(len(sta)):
    n=read(path+sta[k]+'.HNN.mseed')
    e=read(path+sta[k]+'.HNE.mseed')
    z=read(path+sta[k]+'.HNZ.mseed')
    
    #Remove gain
    n[0].data=n[0].data/gain
    e[0].data=e[0].data/gain
    z[0].data=z[0].data/gain
    
    #conver to g
    n[0].data=n[0].data/9.81
    e[0].data=e[0].data/9.81
    z[0].data=z[0].data/9.81
        
    
    pga=r_[abs(n[0].data),abs(e[0].data),abs(z[0].data)].max()
    line='%6s\t%.4f\t%.4f\t%.2f\n' % (sta[k],lonlat[k,0],lonlat[k,1],pga)
    f.write(line)
    
f.close()