from numpy import genfromtxt,zeros,c_,savetxt
from obspy.core import UTCDateTime

f=genfromtxt('/Users/dmelgar/Amatrice2016/afters/afters.txt',delimiter='|')
times=genfromtxt('/Users/dmelgar/Amatrice2016/afters/afters.txt',delimiter='|',usecols=1,dtype='S')
fout='/Users/dmelgar/Amatrice2016/afters/afters.xyz'

seconds=zeros(len(times))
for k in range(len(times)):
    seconds[k]=UTCDateTime(times[k])-UTCDateTime(times[0])
    
days=seconds/86400

#out=c_[f[:,3],f[:,2],days,f[:,10]]
out=c_[f[:,3],f[:,2],f[:,4],f[:,10]]
savetxt(fout,out,fmt='%.4f')