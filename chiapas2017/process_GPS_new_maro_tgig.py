'''
D.Melgar
01/2015
Pre-process Coquimbo GPS data
'''



from glob import glob
from numpy import genfromtxt,where,array,float64,diff,mean,sqrt,sin,cos,arctan2,rad2deg,deg2rad,r_,expand_dims
from obspy import Stream,Trace,read
from os import remove
from obspy.core import UTCDateTime
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees
from datetime import timedelta
from mudpy.forward import lowpass
from matplotlib import pyplot as plt
from shutil import copy
from scipy.signal import filtfilt,butter

path='/Users/dmelgar/Chiapas2017/GPS/new/'
station_file=path+'station_info/gps1Hz_new.sta'
gps_files=glob(path+'raw/kin*')
#gps_files=glob(path+'raw/5Hz/kin*')
time_epi=UTCDateTime('2017-09-08T04:49:21')
epicenter=array([-94.11,14.85,58])  
tmin=timedelta(seconds=10) #Keep these many seconds before p-wave
tmax=timedelta(seconds=180) #Keep these many seconds after p-wave
fcorner=1./10 #Low pass filter corner
dt=0.2 #Sampling itnerval
inter=4.0 #in Hz

make_station_file=False
gps2sac=False
cut_filter=False #Check PGD path
make_plots=True
move_mudpy=False
make_pgd=False



#Some functions to be sued later
def stdecimate(st,factor,order):
    #Anti-alias filter
    b, a = butter(order, 1./factor)
    y = filtfilt(b, a, st[0].data)
    stout=st.copy()
    stout[0].data=y
    #Decimate
    stout[0].decimate(factor,no_filter=True)
    return stout


class WGS84:
    
    #Semimajor axis length (m)
    a = 6378137.0
    #Semiminor axis length (m)
    b = 6356752.3142
    #Ellipsoid flatness (unitless)
    f = (a - b) / a
    #Eccentricity (unitless)
    e = sqrt(f * (2 - f))
    #Speed of light (m/s)
    c = 299792458.
    #Relativistic constant
    F = -4.442807633e-10
    #Earth's universal gravitational constant
    mu = 3.986005e14
    #Earth rotation rate (rad/s)
    omega_ie = 7.2921151467e-5

    def g0(self, L):
        """acceleration due to gravity at the elipsoid surface at latitude L"""
        return 9.7803267715 * (1 + 0.001931851353 * sin(L)**2) / \
                        sqrt(1 - 0.0066943800229 * sin(L)**2)


def ecef2lla(ecef, tolerance=1e-9):
    """Convert Earth-centered, Earth-fixed coordinates to lat, lon, alt.
    Input: ecef - (x, y, z) in (m, m, m)
    Output: lla - (lat, lon, alt) in (decimal degrees, decimal degrees, m)
    """
    from math import atan2,atan
    #Decompose the input
    x = ecef[0]
    y = ecef[1]
    z = ecef[2]
    #Calculate lon
    lon = atan2(y, x)
    #Initialize the variables to calculate lat and alt
    alt = 0
    ellipsoid=WGS84
    N = ellipsoid.a
    p = sqrt(x**2 + y**2)
    lat = 0
    previousLat = 90
    #Iterate until tolerance is reached
    while abs(lat - previousLat) >= tolerance:
        previousLat = lat
        sinLat = z / (N * (1 - ellipsoid.e**2) + alt)
        lat = atan((z + ellipsoid.e**2 * N * sinLat) / p)
        N = ellipsoid.a / sqrt(1 - (ellipsoid.e * sinLat)**2)
        alt = p / cos(lat) - N
    #Return the lla coordinates
    return (rad2deg(lat), rad2deg(lon), alt)


def rotate2neu(x,y,z,lon,lat):
    '''
    Rotate from ECEF 2 local NEU
    '''
    lat=deg2rad(lat)
    lon=deg2rad(lon)
    R=array([[-sin(lat)*cos(lon),-sin(lon)*sin(lat),cos(lat)],[-sin(lon),cos(lon),0],[cos(lon)*cos(lat),cos(lat)*sin(lon),sin(lat)]])
    D=r_[expand_dims(x,0),expand_dims(y,0),expand_dims(z,0)]
    C=R.dot(D)
    n=C[0,:]
    e=C[1,:]
    u=C[2,:]
    
    return n,e,u


if make_station_file:
    fout=open(station_file,'w')
    fout.write('# station,lon,lat,alt(m)\n')
    for k in range(len(gps_files)):
        print k
        f=open(gps_files[k])
        line=f.readline()
        station=gps_files[k].split('_')[-1]
        print station
        f.close()
        coords=genfromtxt(gps_files[k],usecols=[2,3,4],skip_header=4)
        ecef=coords[0,:]
        lat,lon,alt=ecef2lla(ecef)
        fout.write('%s\t%.6f\t%.6f\t%.2f\n' % (station,lon,lat,alt))
    fout.close()
        
        


if gps2sac:
    stanames=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/new/station_info/gps1Hz_new.sta',usecols=0,dtype='S')
    coords=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/new/station_info/gps1Hz_new.sta',usecols=[1,2])
    for k in range(len(gps_files)):

        print k

        #Now read the data
        gps=genfromtxt(gps_files[k],skip_header=4)
        gps=gps[::-1]

        #Initalize obspy stream object
        n=Stream(Trace())
        e=Stream(Trace())
        u=Stream(Trace())
        
        #Fill gaps with zeros
        t=gps[:,1]
        dt=t[1]-t[0]
        print 'dt='+str(dt)
        gap_positions=where(diff(t)>dt+dt/10)[0]+1
        print str(len(gap_positions)+1)+' segments ('+str(len(gap_positions))+' gaps) found'
        if len(gap_positions)>0:  #There are gaps
            for i in range(len(gap_positions)):
                if i==0:        
                    #Fill with data (first trace)
                    n[0].data=gps[0:gap_positions[0],3]
                    e[0].data=gps[0:gap_positions[0],2]
                    u[0].data=gps[0:gap_positions[0],4]
                    #What is first epoch?
                    time=UTCDateTime('2017-09-08T00:00:00')+timedelta(seconds=t[0])-timedelta(seconds=18) #Apply GPS->UTC leapseconds
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
                    n[i].data=gps[gap_positions[i-1]:gap_positions[i],3]
                    e[i].data=gps[gap_positions[i-1]:gap_positions[i],2]
                    u[i].data=gps[gap_positions[i-1]:gap_positions[i],4]
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
            n[i+1].data=gps[gap_positions[i]:,3]
            e[i+1].data=gps[gap_positions[i]:,2]
            u[i+1].data=gps[gap_positions[i]:,4]
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
            n[0].data=gps[:,3]
            e[0].data=gps[:,2]
            u[0].data=gps[:,4]
            #What is first epoch?
            time=UTCDateTime('2017-09-08T00:00:00')+timedelta(seconds=t[0])-timedelta(seconds=18) #Apply GPS->UTC leapseconds
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

        #Write to file
        sta=gps_files[k].split('/')[-1].split('_')[-1]
        n[0].stats.station=sta
        e[0].stats.station=sta
        u[0].stats.station=sta
        n.write(path+'xyz/'+sta+'.LXY.sac',format='SAC')
        e.write(path+'xyz/'+sta+'.LXX.sac',format='SAC')
        u.write(path+'xyz/'+sta+'.LXZ.sac',format='SAC')
        
        #Get coordinates of for that site
        i=where(stanames==sta)[0][0]
        lon=coords[i,0]
        lat=coords[i,1]
        n0,e0,u0=rotate2neu(e[0].data,n[0].data,u[0].data,lon,lat)
        n[0].data=n0
        e[0].data=e0
        u[0].data=u0
        
        n.write(path+'neu/'+sta+'.LXN.sac',format='SAC')
        e.write(path+'neu/'+sta+'.LXE.sac',format='SAC')
        u.write(path+'neu/'+sta+'.LXZ.sac',format='SAC')
        

        
if cut_filter:
    stanames=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/new/station_info/gps1Hz_new.sta',usecols=0,dtype='S')
    coords=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/new/station_info/gps1Hz_new.sta',usecols=[1,2])
    coords[:,0]=coords[:,0]
    for k in range(len(stanames)):
        try:
            sta=stanames[k]
            print sta
            print k
            n=read(path+'neu/'+sta+'.LXN.sac')
            e=read(path+'neu/'+sta+'.LXE.sac')
            u=read(path+'neu/'+sta+'.LXZ.sac')
            
            #Low pass filter
            #n[0].data=lowpass(n[0].data,fcorner,1./n[0].stats.delta,10)
            #e[0].data=lowpass(e[0].data,fcorner,1./e[0].stats.delta,10)
            #u[0].data=lowpass(u[0].data,fcorner,1./u[0].stats.delta,10)
            #decimate
            #n=stdecimate(n,5,4)
            #e=stdecimate(e,5,4)
            #u=stdecimate(u,5,4)
            
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
            ##Interpolate
            n[0].interpolate(sampling_rate=inter)
            e[0].interpolate(sampling_rate=inter)
            u[0].interpolate(sampling_rate=inter)
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
    stanames=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/new/station_info/gps1Hz_new.sta',usecols=0,dtype='S')
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
            plt.savefig(path+'plots/'+sta+'.png')
        except:
            print 'Error on '+sta
        
if move_mudpy:
    stanames=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/station_info/gps5Hz.sta',usecols=0,dtype='S')
    for k in range(len(stanames)):
        try:
            #copy(path+'trim1Hz/'+stanames[k]+'.LXN.sac',path+'trim1Hz/'+stanames[k]+'.disp.n')
            #copy(path+'trim1Hz/'+stanames[k]+'.LXE.sac',path+'trim1Hz/'+stanames[k]+'.disp.e')
            #copy(path+'trim1Hz/'+stanames[k]+'.LXZ.sac',path+'trim1Hz/'+stanames[k]+'.disp.u')
            copy(path+'trim5Hz/'+stanames[k]+'.LXN.sac',path+'trim5Hz/'+stanames[k]+'.disp.n')
            copy(path+'trim5Hz/'+stanames[k]+'.LXE.sac',path+'trim5Hz/'+stanames[k]+'.disp.e')
            copy(path+'trim5Hz/'+stanames[k]+'.LXZ.sac',path+'trim5Hz/'+stanames[k]+'.disp.u')
        except:
            pass
if make_pgd:
    pathout='/Users/dmelgar/PGD/GPS/sac/Melinka2016/'
    tcut=timedelta(minutes=10)
    files=glob(path+'wphase/*LXN.sac')
    for k in range(len(files)):
        sta=files[k].split('/')[-1].split('.')[0]
        n=read(path+'wphase/'+sta+'.LXN.sac')
        e=read(path+'wphase/'+sta+'.LXE.sac')
        u=read(path+'wphase/'+sta+'.LXZ.sac')
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
        
        