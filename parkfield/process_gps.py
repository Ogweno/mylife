from glob import glob
from numpy import genfromtxt,where,array,float64,diff,mean
from obspy import Stream,Trace,read
from os import remove
from obspy.core import UTCDateTime
from datetime import timedelta


path='/Users/dmelgar/Parkfield2004/GPS/'
gps_files=glob(path+'positions/difkin*')
time_epi=UTCDateTime('2004-09-28T17:15:24')
t0=UTCDateTime('2004-09-28T00:00:00')
epicenter=array([-120.374,35.815,7.9]) 
dt=1.0

gps2sac=False
make_pgd=True

if gps2sac:
    stanames=genfromtxt('/Users/dmelgar/Parkfield2004/station_info/parkfield.sta',usecols=0,dtype='S')
    coords=genfromtxt('/Users/dmelgar/Parkfield2004/station_info/parkfield.sta',usecols=[1,2])
    for k in range(len(gps_files)):
        print gps_files[k]
        #Replace AR indicator characters '*' and 'x' with blanks otherwise read breaks down
        f = open(gps_files[k], 'r')
        f1 = open('tmp1', 'w')
        for line in f:
            f1.write(line.replace('*', ''))
        f.close()
        f1.close()
        # 'x'
        f1 = open('tmp1', 'r')
        f2 = open('tmp2', 'w')
        for line in f1:
            f2.write(line.replace('x', ''))
        f1.close()
        f2.close()
        
        #Now read the data
        gps=genfromtxt('tmp2')
        gps=gps[::-1]
        
        #Delete temporary files
        remove('tmp1')
        remove('tmp2')
        
        #Initalize obspy stream object
        n=Stream(Trace())
        e=Stream(Trace())
        u=Stream(Trace())
        
        #Fill gaps with zeros
        t=gps[:,0]
        print 'dt='+str(dt)
        gap_positions=where(diff(t)>dt)[0]+1
        print str(len(gap_positions)+1)+' segments ('+str(len(gap_positions))+' gaps) found'
        if len(gap_positions)>0:  #There are gaps
            for i in range(len(gap_positions)):
                if i==0:        
                    #Fill with data (first trace)
                    n[0].data=gps[0:gap_positions[0],2]/100 #It's in cm, convert to m
                    e[0].data=gps[0:gap_positions[0],1]/100
                    u[0].data=gps[0:gap_positions[0],3]/100
                    #What is first epoch?
                    time=t0+timedelta(seconds=t[0])-timedelta(seconds=13) #Apply GPS->UTC leapseconds
                    #Apply start time
                    n[0].stats.starttime=time
                    e[0].stats.starttime=time
                    u[0].stats.starttime=time
                    #Sampling rate
                    n[0].stats.delta=dt
                    e[0].stats.delta=dt
                    u[0].stats.delta=dt
                else: #it's the second or more trace
                    #Add new trace
                    n+=Trace()
                    e+=Trace()
                    u+=Trace()
                    #Fill with data
                    n[i].data=gps[gap_positions[i-1]:gap_positions[i],2]/100 #It's in cm, convert to m
                    e[i].data=gps[gap_positions[i-1]:gap_positions[i],1]/100
                    u[i].data=gps[gap_positions[i-1]:gap_positions[i],3]/100
                    #Apply start time
                    n[i].stats.starttime=time+timedelta(seconds=t[gap_positions[i-1]])
                    e[i].stats.starttime=time+timedelta(seconds=t[gap_positions[i-1]])
                    u[i].stats.starttime=time+timedelta(seconds=t[gap_positions[i-1]])
                    #Sampling rate
                    n[i].stats.delta=dt
                    e[i].stats.delta=dt
                    u[i].stats.delta=dt
            #Now do last gap to the end of the records
            #Add new trace
            n+=Trace()
            e+=Trace()
            u+=Trace()
            #Fill with data
            n[i+1].data=gps[gap_positions[i]:,2]/100 #It's in cm, convert to m
            e[i+1].data=gps[gap_positions[i]:,1]/100
            u[i+1].data=gps[gap_positions[i]:,3]/100
            #Apply start time
            n[i+1].stats.starttime=time+timedelta(seconds=t[gap_positions[i]])
            e[i+1].stats.starttime=time+timedelta(seconds=t[gap_positions[i]])
            u[i+1].stats.starttime=time+timedelta(seconds=t[gap_positions[i]])
            #Sampling rate
            n[i].stats.delta=dt
            e[i].stats.delta=dt
            u[i].stats.delta=dt
            #Now merge
            n.merge(fill_value='latest')
            e.merge(fill_value='latest')
            u.merge(fill_value='latest')
        else: #No gaps
            n[0].data=gps[:,2]/100 #It's in cm, convert to m
            e[0].data=gps[:,1]/100
            u[0].data=gps[:,3]/100
            #What is first epoch?
            time=t0+timedelta(seconds=t[0])-timedelta(seconds=13) #Apply GPS->UTC leapseconds
            #Apply start time
            n[0].stats.starttime=time
            e[0].stats.starttime=time
            u[0].stats.starttime=time
            #Sampling rate
            n[0].stats.delta=dt
            e[0].stats.delta=dt
            u[0].stats.delta=dt

        #Write to file (twice)
        sta=gps_files[k].split('/')[-1].split('_')[1]
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
        i=where(stanames==sta.upper())[0]
        n[0].stats['sac']['stlo']=coords[i,0]
        n[0].stats['sac']['stla']=coords[i,1]
        e[0].stats['sac']['stlo']=coords[i,0]
        e[0].stats['sac']['stla']=coords[i,1]
        u[0].stats['sac']['stlo']=coords[i,0]
        u[0].stats['sac']['stla']=coords[i,1]
        #Final write
        n.write(path+'proc/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'proc/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'proc/'+sta+'.LXZ.sac',format='SAC')
        
if make_pgd:
    pathout='/Users/dmelgar/PGD/GPS/sac/Parkfield2004/'
    tcut=timedelta(minutes=10)
    files=glob(path+'proc/*LXN.sac')
    for k in range(len(files)):
        sta=files[k].split('/')[-1].split('.')[0]
        n=read(path+'proc/'+sta+'.LXN.sac')
        e=read(path+'proc/'+sta+'.LXE.sac')
        u=read(path+'proc/'+sta+'.LXZ.sac')
        #Trim
        n[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        e[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        u[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        n[0].data=n[0].data-mean(n[0].data[0:11])
        e[0].data=e[0].data-mean(e[0].data[0:11])
        u[0].data=u[0].data-mean(u[0].data[0:11])
        #Write to file
        n.write(pathout+sta+'.LXN.sac',format='SAC')
        e.write(pathout+sta+'.LXE.sac',format='SAC')
        u.write(pathout+sta+'.LXZ.sac',format='SAC')