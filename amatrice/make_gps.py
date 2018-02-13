from numpy import genfromtxt,sqrt

out='/Users/dmelgar/Amatrice2016/GPS/neu_Oct26th/'
stafile='/Users/dmelgar/Amatrice2016/GPS/gps_Oct26th.sta'
threshold=0.005 #In m

#sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/offsets.txt',usecols=0,dtype='S')
#lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/offsets.txt',usecols=[1,2])
#e=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/offsets.txt',usecols=4)
#n=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/offsets.txt',usecols=6)
#u=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/offsets.txt',usecols=8)

#sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
#lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=[1,2])
#e=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=4)
#n=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=6)
#u=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=8)

sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=[1,2])
e=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=4)
n=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=6)
u=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=8)

#Make station file
f=open(stafile,'w')
for k in range(len(sta)):
    line='%s\t%.4f\t%.4f\n' %(sta[k],lonlat[k,0],lonlat[k,1])
    f.write(line)
f.close()

#Make neu files
for k in range(len(sta)):
    
    offset=sqrt((n[k]/1000.)**2+(e[k]/1000.)**2)
    if offset>=threshold:
        f=open(out+sta[k]+'.neu','w')
        f.write('%.6f\n' % (n[k]/1000.))
        f.write('%.6f\n' % (e[k]/1000.))
        f.write('%.6f' % (u[k]/1000.))
        f.close()
