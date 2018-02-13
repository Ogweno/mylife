'''
D.Melgar
01/2015
Pre-process Iquique GPS data
'''



from glob import glob
from numpy import genfromtxt,where,array,float64,diff,mean
from obspy import Stream,Trace,read
from os import remove
from obspy.core import UTCDateTime
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees
from datetime import timedelta
from mudpy.forward import lowpass
from matplotlib import pyplot as plt
from shutil import copy

path='/Users/dmelgar/Amatrice2016/GPS/high_rate/'
gps_files=glob(path+'enu/enu*')
time_epi=UTCDateTime('2016-09-16T01:36:31.50')
epicenter=array([13.2174,42.7123,7.826])  
tmin=timedelta(seconds=10) #Keep these many seconds before p-wave
tmax=timedelta(seconds=30) #Keep these many seconds after p-wave
fcorner=1./10 #Low pass filter corner
dt=1.0 #Sampling itnerval

gps2sac=False
cut_filter=False
make_plots=True
move_mudpy=False
make_pgd=False


if gps2sac:
    stanames=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/stations.txt',usecols=0,dtype='S')
    coords=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/stations.txt',usecols=[1,2])
    for k in range(len(gps_files)):
        print gps_files[k]

        gps=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/enu_'+stanames[k],skip_header=1)
        
        #Fill gaps with zeros
        t=gps[:,0]
        
        n=Stream(Trace())
        e=Stream(Trace())
        u=Stream(Trace())
        
        n[0].data=gps[:,2]/100 #It's in cm, convert to m
        e[0].data=gps[:,1]/100
        u[0].data=gps[:,3]/100
        #What is first epoch?
        time=UTCDateTime('2016-09-16T00:00:00')+timedelta(seconds=t[0])-timedelta(seconds=17) #Apply GPS->UTC leapseconds
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
        i=where(stanames==sta)[0]
        n[0].stats['sac']['stlo']=coords[i,0]
        n[0].stats['sac']['stla']=coords[i,1]
        e[0].stats['sac']['stlo']=coords[i,0]
        e[0].stats['sac']['stla']=coords[i,1]
        u[0].stats['sac']['stlo']=coords[i,0]
        u[0].stats['sac']['stla']=coords[i,1]
        #Final write
        n[0].write(path+'proc/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'proc/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'proc/'+sta+'.LXZ.sac',format='SAC')
        
if cut_filter:
    stanames=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/stations.txt',usecols=0,dtype='S')
    coords=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/stations.txt',usecols=[1,2])
    for k in range(len(stanames)):
        try:
            sta=stanames[k]
            print sta
            print k
            n=read(path+'proc/'+sta+'.LXN.sac')
            e=read(path+'proc/'+sta+'.LXE.sac')
            u=read(path+'proc/'+sta+'.LXZ.sac')
            #Low pass filter
            #n[0].data=lowpass(n[0].data,fcorner,1./n[0].stats.delta,10)
            #e[0].data=lowpass(e[0].data,fcorner,1./e[0].stats.delta,10)
            #u[0].data=lowpass(u[0].data,fcorner,1./u[0].stats.delta,10)
            #Get station to hypocenter delta distance
            delta=locations2degrees(coords[k,0],coords[k,1],epicenter[1],epicenter[0])
            #Get p-time to site
            tt=getTravelTimes(delta,epicenter[2])
            tp=timedelta(seconds=float64(tt[0]['time']))
            #Trim
            n[0].trim(starttime=time_epi+tp-tmin,endtime=time_epi+tp+tmax)
            e[0].trim(starttime=time_epi+tp-tmin,endtime=time_epi+tp+tmax)
            u[0].trim(starttime=time_epi+tp-tmin,endtime=time_epi+tp+tmax)
            #Remove first epoch
            n[0].data=n[0].data-n[0].data[0]
            e[0].data=e[0].data-e[0].data[0]
            u[0].data=u[0].data-u[0].data[0]
            #Write to file
            n.write(path+'trim/'+sta+'.LXN.sac',format='SAC')
            e.write(path+'trim/'+sta+'.LXE.sac',format='SAC')
            u.write(path+'trim/'+sta+'.LXZ.sac',format='SAC')
        except:
            print 'Error on '+sta
        
if make_plots:
    stanames=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/high_rate/enu/stations.txt',usecols=0,dtype='S')
    for k in range(len(stanames)):
        print stanames[k]
        sta=stanames[k]
        n=read(path+'trim/'+sta+'.LXN.sac')
        e=read(path+'trim/'+sta+'.LXE.sac')
        z=read(path+'trim/'+sta+'.LXZ.sac')
        plt.figure()
        plt.subplot(311)
        plt.plot(n[0].times(),n[0].data,'r')
        plt.grid()
        plt.ylabel('North (m)')
        plt.title('Station '+stanames[k])
        plt.subplot(312)
        plt.plot(e[0].times(),e[0].data,'b')
        plt.grid()
        plt.ylabel('East (m)')
        plt.subplot(313)
        plt.plot(z[0].times(),z[0].data,'g')
        plt.grid()
        plt.ylabel('Up (m)')
        plt.savefig(path+'plots/'+sta+'.png')
        
if move_mudpy:
    stanames=genfromtxt('/Users/dmelgar/Slip_inv/iquique_sm/data/station_info/gps.gflist',usecols=0,dtype='S')
    for k in range(len(stanames)):
        move(path+'filt/'+stanames[k]+'.LXN.sac','/Users/dmelgar/Slip_inv/iquique_gps/data/waveforms/'+stanames[k]+'.disp.n')
        move(path+'filt/'+stanames[k]+'.LXE.sac','/Users/dmelgar/Slip_inv/iquique_gps/data/waveforms/'+stanames[k]+'.disp.e')
        move(path+'filt/'+stanames[k]+'.LXZ.sac','/Users/dmelgar/Slip_inv/iquique_gps/data/waveforms/'+stanames[k]+'.disp.u')
        
if make_pgd:
    pathout='/Users/dmelgar/PGD/GPS/sac/Iquique2014/'
    time_epi=UTCDateTime('2014-04-01T23:46:47')
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
        
        