from numpy import genfromtxt,where,c_,savetxt,zeros

f=genfromtxt('/Users/dmelgar/DART_analysis/stations_degminsec.txt',usecols=[0,1,2,3,4,5,6,7,8,9])
ns=genfromtxt('/Users/dmelgar/DART_analysis/stations_degminsec.txt',usecols=[5],dtype='S')
ew=genfromtxt('/Users/dmelgar/DART_analysis/stations_degminsec.txt',usecols=[9],dtype='S')
sta=genfromtxt('/Users/dmelgar/DART_analysis/stations_degminsec.txt',usecols=[0],dtype='f')
fout='/Users/dmelgar/DART_analysis/stations.txt'
fout2='/Users/dmelgar/DART_analysis/stations_filtered.txt'


lat=f[:,2]+f[:,3]/60.+f[:,4]/3600.
lon=f[:,6]+f[:,7]/60.+f[:,8]/3600.

i=where(ns=='S')[0]
lat[i]=-lat[i]

i=where(ew=='W')[0]
lon[i]=-lon[i]

out=c_[sta,lon,lat]
savetxt(fout,out,fmt='%d\t%.4f\t%.4f')


#filter
keep=genfromtxt('/Users/dmelgar/DART_analysis/keep.txt')
out2=zeros((len(keep),3))
for k in range(len(keep)):
    i=where(sta==keep[k])[0]
    out2[k,:]=out[i,:]
    
savetxt(fout2,out2,fmt='%d\t%.4f\t%.4f')
