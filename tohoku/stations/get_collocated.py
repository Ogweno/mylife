from numpy import array,r_,where
from obspy.core.util.geodetics import gps2DistAzimuth

pbo_list='/Users/dmelgar/KMLs/pbo_existing.txt'
sm_list='/Users/dmelgar/KMLs/new_sm.txt'
bb_sm_list='/Users/dmelgar/KMLs/new_bb_sm.txt'
dist_filter=5000
fout1='/Users/dmelgar/KMLs/pbo_new_sm_collocated'+str(dist_filter)+'.txt'
fout2='/Users/dmelgar/KMLs/pbo_new_bb_sm_collocated'+str(dist_filter)+'.txt'

#Parse PBO
f=open(pbo_list,'r')
pbo_sta=[]
pbo_lon=array([])
pbo_lat=array([])
while True:
    line=f.readline()
    if '>' in line:
        pbo_sta.append(line.split('"')[1])
    else:
        try:
            pbo_lon=r_[pbo_lon,float(line.split()[0])]
            pbo_lat=r_[pbo_lat,float(line.split()[1])]
        except:
            pass
    if not line:
        break
        
#Parse Proposed stations
f=open(sm_list,'r')
sm_sta=[]
sm_lon=array([])
sm_lat=array([])
while True:
    line=f.readline()
    if '>' in line:
        sm_sta.append(line.split('"')[1])
    else:
        try:
            sm_lon=r_[sm_lon,float(line.split()[0])]
            sm_lat=r_[sm_lat,float(line.split()[1])]
        except:
            pass
    if not line:
        break   
        
f=open(bb_sm_list,'r')
bb_sm_sta=[]
bb_sm_lon=array([])
bb_sm_lat=array([])
while True:
    line=f.readline()
    if '>' in line:
        bb_sm_sta.append(line.split('"')[1])
    else:
        try:
            bb_sm_lon=r_[bb_sm_lon,float(line.split()[0])]
            bb_sm_lat=r_[bb_sm_lat,float(line.split()[1])]
        except:
            pass
    if not line:
        break   
        
#Compute distances
i=where((pbo_lat<42) & (pbo_lat>35))[0]
pbo_lon=pbo_lon[i]
pbo_lat=pbo_lat[i]

f1=open(fout1,'w')
f2=open(fout2,'w')
for kpbo in range(len(pbo_lon)):
    print kpbo
    for kseis in range(len(sm_lon)):
        d,az,baz=gps2DistAzimuth(pbo_lat[kpbo],pbo_lon[kpbo],sm_lat[kseis],sm_lon[kseis])
        if d< dist_filter:
            print 'Found collocation! '+pbo_sta[i[kpbo]]+' and '+sm_sta[kseis]+' are '+str(d)+'m apart'
            out='%s\t%12.6f\t%12.6f\t%s\t%12.6f\t%12.6f\t%7.1f\n' %(pbo_sta[i[kpbo]],pbo_lon[kpbo],pbo_lat[kpbo],sm_sta[kseis],sm_lon[kseis],sm_lat[kseis],d)
            f1.write(out)
    for kseis in range(len(bb_sm_lon)):
        d,az,baz=gps2DistAzimuth(pbo_lat[kpbo],pbo_lon[kpbo],bb_sm_lat[kseis],bb_sm_lon[kseis])
        if d< dist_filter:
            print 'Found collocation! '+pbo_sta[i[kpbo]]+' and '+bb_sm_sta[kseis]+' are '+str(d)+'m apart'
            out='%s\t%12.6f\t%12.6f\t%s\t%12.6f\t%12.6f\t%7.1f\n' %(pbo_sta[i[kpbo]],pbo_lon[kpbo],pbo_lat[kpbo],bb_sm_sta[kseis],bb_sm_lon[kseis],sm_lat[kseis],d)
            f2.write(out)
f1.close()
f2.close()