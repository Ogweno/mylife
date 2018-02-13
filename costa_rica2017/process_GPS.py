'''
D.Melgar

Pre-process GPS data
'''



from glob import glob
from numpy import genfromtxt,where,array,float64,diff,mean,sqrt,sin,cos,arctan2,rad2deg,deg2rad,r_,expand_dims
from obspy import Stream,Trace,read
from obspy.core import UTCDateTime
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees
from datetime import timedelta
from mudpy.forward import lowpass
from matplotlib import pyplot as plt
from shutil import copy
from scipy.signal import filtfilt,butter

path='/Users/dmelgar/Costa_Rica2017/GPS/'
station_file=path+'station_info/gps.sta'
gps_files=glob(path+'raw/enu*')
time_epi=UTCDateTime('2017-11-13T02:28:23')
epicenter=array([-84.487,9.515,19.4])  
tmin=timedelta(seconds=10) #Keep these many seconds before p-wave
tmax=timedelta(seconds=90) #Keep these many seconds after p-wave
leap_seconds=18
event_day='2017-11-13T00:00:00'

# What do you want to do?
gps2sac=False
cut=False
make_plots=True


#Convert from text to SAC
if gps2sac:
    stanames=genfromtxt(station_file,usecols=0,dtype='S')
    coords=genfromtxt(station_file,usecols=[1,2])
    for k in range(len(gps_files)):

        print k

        #Now read the data
        gps=genfromtxt(gps_files[k],skip_header=4)

        #Initalize obspy stream object
        n=Stream(Trace())
        e=Stream(Trace())
        u=Stream(Trace())
        
        #Fill gaps with zeros
        t=gps[:,0]
        dt=t[1]-t[0]
        print 'dt='+str(dt)
        gap_positions=where(diff(t)>dt+dt/10)[0]+1
        print str(len(gap_positions)+1)+' segments ('+str(len(gap_positions))+' gaps) found'
        if len(gap_positions)>0:  #There are gaps
            for i in range(len(gap_positions)):
                if i==0:        
                    #Fill with data (first trace)
                    n[0].data=gps[0:gap_positions[0],2]
                    e[0].data=gps[0:gap_positions[0],1]
                    u[0].data=gps[0:gap_positions[0],3]
                    #What is first epoch?
                    time=UTCDateTime(event_day)+timedelta(seconds=t[0])-timedelta(seconds=leap_seconds) #Apply GPS->UTC leapseconds
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
                    n[i].data=gps[gap_positions[i-1]:gap_positions[i],2]
                    e[i].data=gps[gap_positions[i-1]:gap_positions[i],1]
                    u[i].data=gps[gap_positions[i-1]:gap_positions[i],3]
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
            n[i+1].data=gps[gap_positions[i]:,2]
            e[i+1].data=gps[gap_positions[i]:,1]
            u[i+1].data=gps[gap_positions[i]:,3]
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
            n[0].data=gps[:,2]
            e[0].data=gps[:,1]
            u[0].data=gps[:,3]
            #What is first epoch?
            time=UTCDateTime(event_day)+timedelta(seconds=t[0])-timedelta(seconds=leap_seconds) #Apply GPS->UTC leapseconds
            #Apply start time
            n[0].stats.starttime=time
            e[0].stats.starttime=time
            u[0].stats.starttime=time
            #Sampling rate
            n[0].stats.delta=dt
            e[0].stats.delta=dt
            u[0].stats.delta=dt
        #remove first sample
        n[0].data=n[0].data-n[0].data[0]
        e[0].data=e[0].data-e[0].data[0]
        u[0].data=u[0].data-u[0].data[0]
        
        #convert from cm -> m
        n[0].data=n[0].data/100.
        e[0].data=e[0].data/100.
        u[0].data=u[0].data/100.

        #Write to file
        sta=gps_files[k][-4:]
        n[0].stats.station=sta
        e[0].stats.station=sta
        u[0].stats.station=sta
        n.write(path+'sac/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'sac/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'sac/'+sta+'.LXZ.sac',format='SAC')
        
        

 #Trim around time of earthquake       
if cut:
    stanames=genfromtxt(station_file,usecols=0,dtype='S')
    coords=genfromtxt(station_file,usecols=[1,2])
    for k in range(len(stanames)):
        try:
            sta=stanames[k]
            print sta
            print k
            n=read(path+'sac/'+sta+'.LXN.sac')
            e=read(path+'sac/'+sta+'.LXE.sac')
            u=read(path+'sac/'+sta+'.LXZ.sac')
            
            #Get station to hypocenter delta distance
            delta=locations2degrees(coords[k,1],coords[k,0],epicenter[1],epicenter[0])
            
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

            #Add headers
            n[0].stats['sac']['stlo']=coords[k,0]
            n[0].stats['sac']['stla']=coords[k,1]
            e[0].stats['sac']['stlo']=coords[k,0]
            e[0].stats['sac']['stla']=coords[k,1]
            u[0].stats['sac']['stlo']=coords[k,0]
            u[0].stats['sac']['stla']=coords[k,1]
            n[0].stats['sac']['evlo']=epicenter[0]
            n[0].stats['sac']['evla']=epicenter[1]
            e[0].stats['sac']['evlo']=epicenter[0]
            e[0].stats['sac']['evla']=epicenter[1]
            u[0].stats['sac']['evlo']=epicenter[0]
            u[0].stats['sac']['evla']=epicenter[1]       

            #Write to file
            n.write(path+'trim/'+sta+'.LXN.sac',format='SAC')
            e.write(path+'trim/'+sta+'.LXE.sac',format='SAC')
            u.write(path+'trim/'+sta+'.LXZ.sac',format='SAC')

        except:
            print 'Error on '+sta
        
if make_plots:
    stanames=genfromtxt(station_file,usecols=0,dtype='S')
    for k in range(len(stanames)):
        print stanames[k]
        try:
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
            plt.savefig(path+'plots/short/'+sta+'.png')
        except:
            print 'Error on '+sta
        
