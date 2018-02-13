from numpy import genfromtxt,where


sta_picks=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=0,dtype='S')
pick_type=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=2,dtype='S')
pick_weight=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=3)
hour=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=4)
minute=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=5)
second=genfromtxt('/Users/dmelgar/Chiapas2017/locate/picks.txt',usecols=6)

ssn_sta=genfromtxt(u'/Users/dmelgar/SSN/all_stations.txt',usecols=1,dtype='S')
ssn_lon=genfromtxt(u'/Users/dmelgar/SSN/all_stations.txt',usecols=2)
ssn_lat=genfromtxt(u'/Users/dmelgar/SSN/all_stations.txt',usecols=3)

f=open('/Users/dmelgar/Chiapas2017/locate/elocate.dat','w')

count=1
for k in range(len(sta_picks)):
    i=where(ssn_sta==sta_picks[k])[0]
    lon=ssn_lon[i]
    lat=ssn_lat[i]
    print sta_picks[k]+' '+str(i)
    
    if 'P' in pick_type[k]:
        pt='P'
        chan='Z'
    else:
        pt='S'
        chan='E'
        
    if 'I' in pick_type[k]:
        ei='i'
    else:
        ei='e'
    
    
    if len(i)>0:
        # LSS    Z  2016 08 24 01 36 36.831        1 i P         0 -         1   42.55820  12.9689      0.          0  9999999
        line='%4s   %s  2017 09 19 18 %d %6.3f       %2d %s %s         0 -         %2d  %8.5f %8.5f      0.          0  9999999\n' % (sta_picks[k],chan,minute[k],second[k],count,ei,pt,count,lat,lon)
        f.write(line)
        count+=1
f.close()