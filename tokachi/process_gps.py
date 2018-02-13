from glob import glob
from numpy import genfromtxt,where,array,mean
from obspy import Stream,Trace,read
from obspy.core import UTCDateTime
from datetime import timedelta
from string import rjust

path='/Users/dmelgar/Tokachi2003/GPS/'
gps_files='/Users/dmelgar/Tokachi2003/GPS/positions/hokkneu.txt'
time_epi=UTCDateTime('2003-09-25T19:50:06')
t0=UTCDateTime('2003-09-25T19:50:06')
epicenter=array([143.904,41.775,27]) 
dt=1.0

gps2sac=True
make_pgd=True

if gps2sac:
    stanames=genfromtxt('/Users/dmelgar/Tokachi2003/GPS/station_info/tokachi.sta',usecols=0,dtype='i')
    coords=genfromtxt('/Users/dmelgar/Tokachi2003/GPS/station_info/tokachi.sta',usecols=[1,2])
    pos=genfromtxt(gps_files)
    for k in range(len(stanames)):        
        #Initalize obspy stream object
        n=Stream(Trace())
        e=Stream(Trace())
        u=Stream(Trace())
        #Find data
        ista=where(pos[:,0]==stanames[k])[0]
        t=pos[ista,4]
        n[0].data=pos[ista,1] #It's in cm, convert to m
        e[0].data=pos[ista,2]
        u[0].data=pos[ista,3]
        #What is first epoch?
        time=t0+timedelta(seconds=t[0]) #Apply GPS->UTC leapseconds
        #Apply start time
        n[0].stats.starttime=time
        e[0].stats.starttime=time
        u[0].stats.starttime=time
        #Sampling rate
        n[0].stats.delta=dt
        e[0].stats.delta=dt
        u[0].stats.delta=dt
        #Write to file (twice)
        sta=rjust(str(stanames[k]),4,'0')
        n[0].stats.station=sta
        e[0].stats.station=sta
        u[0].stats.station=sta
        n.write(path+'proc/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'proc/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'proc/'+sta+'.LXZ.sac',format='SAC')
        #Re-read
        n=read(path+'proc/'+sta+'.LXN.sac')
        e=read(path+'proc/'+sta+'.LXE.sac')
        u=read(path+'proc/'+sta+'.LXZ.sac')
        #get positions
        n[0].stats['sac']['stlo']=coords[k,0]
        n[0].stats['sac']['stla']=coords[k,1]
        e[0].stats['sac']['stlo']=coords[k,0]
        e[0].stats['sac']['stla']=coords[k,1]
        u[0].stats['sac']['stlo']=coords[k,0]
        u[0].stats['sac']['stla']=coords[k,1]
        #Final write
        n.write(path+'proc/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'proc/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'proc/'+sta+'.LXZ.sac',format='SAC')
        
if make_pgd:
    pathout='/Users/dmelgar/PGD/GPS/sac/Tokachi2003/'
    tcut=timedelta(minutes=6)
    files=glob(path+'proc/*LXN.sac')
    for k in range(len(files)):
        sta=files[k].split('/')[-1].split('.')[0]
        print sta
        n=read(path+'proc/'+sta+'.LXN.sac')
        e=read(path+'proc/'+sta+'.LXE.sac')
        u=read(path+'proc/'+sta+'.LXZ.sac')
        #Trim
        n[0].trim(starttime=time_epi,endtime=n[0].stats.endtime)
        e[0].trim(starttime=time_epi,endtime=e[0].stats.endtime)
        u[0].trim(starttime=time_epi,endtime=u[0].stats.endtime)
        n[0].data=n[0].data-mean(n[0].data[0:9])
        e[0].data=e[0].data-mean(e[0].data[0:9])
        u[0].data=u[0].data-mean(u[0].data[0:9])
        #Write to file
        n.write(pathout+sta+'.LXN.sac',format='SAC')
        e.write(pathout+sta+'.LXE.sac',format='SAC')
        u.write(pathout+sta+'.LXZ.sac',format='SAC')