from numpy import genfromtxt,where

g1=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat')
sta1=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
g2=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat')
sta2=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
combination_file='/Users/dmelgar/Amatrice2016/GPS/Oct26-30_combined.txt'

f=open(combination_file,'w')
f.write('# sta,lon,lat,east(m),north(m),up(m)\n')
for k in range(len(g1)):
    sta=sta1[k]
    i=where(sta2==sta)[0]
    if len(i)!=0:
        lon=g1[k,1]
        lat=g1[k,2]
        east=(g1[k,4]+g2[i,4])/1000
        north=(g1[k,6]+g2[i,6])/1000
        up=(g1[k,8]+g2[i,8])/1000
        f.write('%s\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\n' %(sta,lon,lat,east,north,up))
f.close()
    