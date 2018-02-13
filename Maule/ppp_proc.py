from glob import glob
from numpy import genfromtxt,r_,savetxt,zeros,fliplr
from gpstools import ecef2lla
from obspy.core import UTCDateTime
from obspy import Stream,Trace
from datetime import timedelta

ppp_path=u'/Users/dmelgar/Maule2010/GPS/PPPAR/'


#GPS day
gps_day=UTCDateTime('2010-02-27T00:00:00')
leapseconds=timedelta(seconds=15)

##Make station list (use first position)
#List files
#f=glob(ppp_path+'pos*')
#out=zeros((len(f),3))
#fout=open(ppp_path+'ppp_list_forkml.txt','w')
#fout.write('# Lon    Lat    height(m)\n')
#for k in range(len(f)):
#    sta=f[k].split('_')[-1]
#    pos=genfromtxt(f[k],usecols=[2,3,4])
#    pos=pos[0,:]
#    lon,lat,z=ecef2lla(pos)
#    lineout='%14.6f\t%14.6f\t%14.6f\t%s\n' %(lon,lat,z,sta)
#    lineout='%14.6f\t%14.6f\t%s\n' %(lon,lat,sta)
#    fout.write(lineout)
#fout.close()

#Conver to SAC
f=glob(ppp_path+'difkin_????')
for k in range(len(f)):
    sta=f[k].split('_')[-1]
    g=genfromtxt(f[k],usecols=[0,1,2,3])
    t=g[:,0]
    e=g[:,1]/100
    n=g[:,2]/100
    u=g[:,3]/100
    eout=Stream(Trace())
    nout=Stream(Trace())
    uout=Stream(Trace())
    eout[0].data=e[::-1]
    nout[0].data=n[::-1]
    uout[0].data=u[::-1]
    #What time is it?
    eout[0].stats.starttime=gps_day+timedelta(seconds=t[0])-leapseconds
    nout[0].stats.starttime=gps_day+timedelta(seconds=t[0])-leapseconds
    uout[0].stats.starttime=gps_day+timedelta(seconds=t[0])-leapseconds
    eout[0].stats.station=sta
    nout[0].stats.station=sta
    uout[0].stats.station=sta
    delta=t[0]-t[1]
    eout[0].stats.delta=delta
    nout[0].stats.delta=delta
    uout[0].stats.delta=delta
    #Write
    eout.write(ppp_path+'SAC/'+sta+'.LXE.sac',format='SAC')
    nout.write(ppp_path+'SAC/'+sta+'.LXN.sac',format='SAC')
    uout.write(ppp_path+'SAC/'+sta+'.LXZ.sac',format='SAC')
    
    
