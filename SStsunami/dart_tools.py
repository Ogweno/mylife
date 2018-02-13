'''
Tools for manipulating DART data
'''

def download_all(station_file,out_path):
    '''
    Download all avialble DART data
    '''
    
    from numpy import genfromtxt,arange
    import urllib
    from os import path,makedirs
    
    stations=genfromtxt(station_file,usecols=0,dtype='i')
    years=arange(2000,2016)
    
    for ksta in range(len(stations)):
        
        print 'Station '+str(stations[ksta])
        folder=out_path+str(stations[ksta])+'/'
        
        #Check if station directory exists, if not then make it
        if path.exists(folder):
            pass
        else:
            makedirs(folder)
        
        #Loop over years
        for kyr in range (len(years)):
            
            print '... year '+str(years[kyr])
            
            #Form url
            url='http://www.ndbc.noaa.gov/view_text_file.php?filename='+str(stations[ksta])+'t'+str(years[kyr])+'.txt.gz&dir=data/historical/dart/'
            out_file=folder+str(stations[ksta])+'.'+str(years[kyr])+'.dart'
            
            #Get file
            urllib.urlretrieve(url,out_file)
            
        #Get march 2016
        url='http://www.ndbc.noaa.gov/view_text_file.php?filename='+str(stations[ksta])+'32016.txt.gz&dir=data/dart/Mar/'
        out_file=folder+str(stations[ksta])+'.2016.dart'
        urllib.urlretrieve(url,out_file)
        
    return 
            

def cleanup_empty(path):
    '''
    Delete empty files
    '''
    
    from glob import glob
    from os import remove
    
    folders=glob(path+'*')
    for k in range(len(folders)):
        files=glob(folders[k]+'/*.dart')
        for kf in range(len(files)):
            f=open(files[kf],'r')
            if 'Unable to access data file' in f.readline():
                f.close()
                print '... deleting '+files[kf]
                remove(files[kf])
            else:
                f.close()      
    return
    

def dart2sac(dart_file,dt_in_sec=15,ioc=False):
    '''
    Make one .dart text file to sac
    '''
    
    from glob import glob
    from numpy import genfromtxt,where,zeros,arange
    import datetime
    from scipy.interpolate import interp1d
    from obspy import Stream,Trace
    from os import path as p
    
    
    #Station?
    sta=dart_file.split('/')[-1].split('.')[0]
    
    dart=genfromtxt(dart_file,skip_header=1)
    #Parse data to individual variables
    yr=dart[:,0].astype(int)
    mo=dart[:,1].astype(int)
    dy=dart[:,2].astype(int)
    hr=dart[:,3].astype(int)
    mi=dart[:,4].astype(int)
    se=dart[:,5].astype(int)
    if ioc==False:
        eta=dart[:,7]
    else:
        eta=dart[:,6]  
    
    #Exclude no data samples (no data=9999)
    i=where(eta!=9999)[0]
    
    #Get rid of 'em
    yr=yr[i] ; mo=mo[i] ; dy=dy[i]
    hr=hr[i] ; mi=mi[i] ; se=se[i]
    eta=eta[i]
    
    ifix=where(hr==24)[0]
    hr[ifix]=0
    
    #Make array of seconds since first epoch in record
    if len(yr)>0:
        t=zeros(len(eta))
        for kt in range(len(yr)):
            #Get first epoch
            if kt==0:
                t1=datetime.datetime(yr[kt],mo[kt],dy[kt],hr[kt],mi[kt],se[kt])
            else:
                tnow=datetime.datetime(yr[kt],mo[kt],dy[kt],hr[kt],mi[kt],se[kt])
                t[kt]=(tnow-t1).total_seconds()
                
        #resample to regular interval
        t_out=arange(0,t[-1],dt_in_sec)
        f=interp1d(t,eta)
        eta_out=f(t_out)
        
        #Put in obspy object and save as SAC
        st=Stream(Trace())
        st[0].data=eta_out
        st[0].stats.delta=dt_in_sec
        st[0].stats.starttime=t1
        st[0].stats.station=sta  
        
    return st
        
          

def make_sac(path,dt_in_sec=15.):
    '''
    Resample to regular interval and Convert all files to sac
    '''
    
    from glob import glob
    from numpy import genfromtxt,where,zeros,arange
    import datetime
    from scipy.interpolate import interp1d
    from obspy import Stream,Trace
    from os import path as p
    
    #Loop over stations
    folders=glob(path+'*')
    for k in range(21,len(folders)):
        
        #Loop over years in that station
        files=glob(folders[k]+'/*.dart')
        base_folder=p.dirname(files[0])+'/'
        for kf in range(len(files)):
            
            print files[kf]
            
            #Station?
            sta=files[kf].split('/')[-1].split('.')[0]
            year=files[kf].split('/')[-1].split('.')[1]
            
            if int(year)>2006:
                dart=genfromtxt(files[kf])
            else:
                dart=genfromtxt(files[kf],skip_header=1)
            #Parse data to individual variables
            yr=dart[:,0].astype(int)
            mo=dart[:,1].astype(int)
            dy=dart[:,2].astype(int)
            hr=dart[:,3].astype(int)
            mi=dart[:,4].astype(int)
            se=dart[:,5].astype(int)
            eta=dart[:,7]
            
            #Exclude no data samples (no data=9999)
            i=where(eta!=9999)[0]
            
            #Get rid of 'em
            yr=yr[i] ; mo=mo[i] ; dy=dy[i]
            hr=hr[i] ; mi=mi[i] ; se=se[i]
            eta=eta[i]
            
            #Make array of seconds since first epoch in record
            if len(yr)>0:
                t=zeros(len(eta))
                for kt in range(len(yr)):
                    #Get first epoch
                    if kt==0:
                        t1=datetime.datetime(yr[kt],mo[kt],dy[kt],hr[kt],mi[kt],se[kt])
                    else:
                        tnow=datetime.datetime(yr[kt],mo[kt],dy[kt],hr[kt],mi[kt],se[kt])
                        t[kt]=(tnow-t1).total_seconds()
                        
                #resample to regular interval
                t_out=arange(0,t[-1],dt_in_sec)
                f=interp1d(t,eta)
                eta_out=f(t_out)
                
                #Put in obspy object and save as SAC
                st=Stream(Trace())
                st[0].data=eta_out
                st[0].stats.delta=dt_in_sec
                st[0].stats.starttime=t1
                st[0].stats.station=sta
                
                #Write
                st.write(base_folder+sta+'.'+year+'.sac',format='SAC')

                        
def match_events_and_stations(catalogue_file,stations_file,station_path,distance=2000.):
    '''
    Analyze CMT events file, look at stations file and see whoch stations are 
    within a certain threshold distance
    '''
    
    from pyproj import Geod
    from numpy import genfromtxt,ones,where,zeros
    from glob import glob
    from datetime import datetime
    
    #Read station coordinates
    stations=genfromtxt(stations_file)
    sta=stations[:,0].astype(int)
    sta_lon=stations[:,1]
    sta_lat=stations[:,2]
    
    #Read events and get coordinates
    cmt_dates=genfromtxt(catalogue_file,usecols=12,dtype='S')
    cmt_lon=genfromtxt(catalogue_file,usecols=0)
    cmt_lat=genfromtxt(catalogue_file,usecols=1)
    
    #Parse dates
    year=zeros(len(cmt_dates))
    month=zeros(len(cmt_dates))
    day=zeros(len(cmt_dates))
    hour=zeros(len(cmt_dates))
    minute=zeros(len(cmt_dates))
    
    for k in range(len(cmt_dates)):
        if len(cmt_dates[k])>10:
            date_string=cmt_dates[k]
            year[k]=int(date_string[0:4])
            month[k]=int(date_string[4:6])
            day[k]=int(date_string[6:8])
            hour[k]=int(date_string[8:10])
            minute[k]=int(date_string[10:12])
    
    #Look only at events after 2005
    i=where(year>=2005)[0]
    year=year[i] ; month=month[i] ; day=day[i]
    hour=hour[i] ;minute=minute[i]
    cmt_lon=cmt_lon[i] ; cmt_lat=cmt_lat[i]
    
    # Event loop
    g=Geod(ellps='WGS84')
    
    for kevent in range(len(i)):
        
        #replicate single event lat/lon by number of stations
        ev_lon=cmt_lon[kevent]*ones(len(sta_lon))
        ev_lat=cmt_lat[kevent]*ones(len(sta_lon))
        az,baz,dist=g.inv(ev_lon,ev_lat,sta_lon,sta_lat)
        dist=dist/1000
        # How many events within distance?
        Ndist=where(dist<=distance)[0]
        
        #Eevent origin time object
        if kevent==0:
            event_origin=[datetime(int(year[kevent]),int(month[kevent]),int(day[kevent]),int(hour[kevent]),int(minute[kevent]))]
        else:
            event_origin.append(datetime(int(year[kevent]),int(month[kevent]),int(day[kevent]),int(hour[kevent]),int(minute[kevent])))
        
        
        if len(Ndist)==0:
            event_stations=['None']
            if kevent==0:
                distances=['None']
            else:
                distances.append(['None'])
        else:
            event_stations=[]
            station_dist=[]
            for ksta in range(len(Ndist)):
                station_file=glob(station_path+str(sta[Ndist[ksta]])+'/'+str(sta[Ndist[ksta]])+'.'+str(int(year[kevent]))+'.sac')
                if station_file!=[]:
                    event_stations.append(station_file)
                    #And save distance info
                    station_dist.append(dist[Ndist[ksta]])
            if kevent==0:
                distances=station_dist[:]
            else:
                distances.append(station_dist)
                    
            if event_stations==[]:
                event_stations=['None']
        
        if kevent==0:
            stations_with_data=event_stations[:]
        else:
            stations_with_data.append(event_stations)
            
    return stations_with_data,event_origin,distances
    

def station_event_plots(stations_with_data,event_origin,distances,out_path,twindow_in_hrs=4.0,fcorner=1./3600,spike_level=4.0):
    '''
    Plots of stations for each event
    '''
    
    from obspy import read
    from matplotlib import pyplot as plt 
    from obspy.core import UTCDateTime
    from datetime import timedelta
    from mudpy.forward import highpass
    from string import rjust
    from numpy import where
                
    for k in range(len(stations_with_data)):
        
        if stations_with_data[k]=='None' or stations_with_data[k][0]=='None':
            pass
        else:
            files=stations_with_data[k]
            event_start=event_origin[k].strftime('%Y-%m-%dT%H:%M')
            print event_start
            event_start=UTCDateTime(event_start)
            event_end=event_start+timedelta(hours=twindow_in_hrs)
            #Event to station distances
            dist=distances[k]
            
            #init figure
            fig, axarr = plt.subplots(len(files), 1)
            
            for ksta in range(len(files)):
                
                try:
                    len(axarr)
                    ax=axarr[ksta]
                except:
                    ax=axarr
                station=files[ksta][0].split('/')[-1].split('.')[0]
                
                st=read(files[ksta][0])
                st[0].trim(starttime=event_start,endtime=event_end)

                if st[0].stats.npts>0:
                    st[0].data=st[0].data-st[0].data.mean()
                    #Despike
                    i=where(st[0].data>spike_level)[0]
                    st[0].data[i]=0
                    st[0].data=highpass(st[0].data,fcorner,1./st[0].stats.delta,2)
                    ax.plot(st[0].times()/3600.,st[0].data,'k')
                    yl=[min(st[0].data),max(st[0].data)]
                    #Reference lines
                    ax.plot([dist[ksta]/(4.*3600),dist[ksta]/(4.*3600)],[-50,50],'r--')
                    ax.plot([dist[ksta]/(6.*3600),dist[ksta]/(6.*3600)],[-50,50],'r--')
                    ax.plot([dist[ksta]/(8.*3600),dist[ksta]/(8.*3600)],[-50,50],'r--')
                    ax.plot([dist[ksta]/617.,dist[ksta]/617.],[-50,50],'b--')
                    ax.plot([dist[ksta]/797.,dist[ksta]/797.],[-50,50],'b--')
                    ax.plot([dist[ksta]/943.,dist[ksta]/943.],[-50,50],'b--')
                    ax.set_ylim(yl)
                    
                else:
                    ax.plot(0,0,'k')
                ax.legend(['station %s, d=%skm' % (station,int(dist[ksta]))],frameon=False)
                ax.set_ylabel('m')
            if ksta==0:
                ax.set_title('Event at '+str(event_start))
            if ksta==len(files)-1:
                ax.set_xlabel('Hours since OT')
            event_string=str(event_start.year)+rjust(str(event_start.month),2,'0')+rjust(str(event_start.day),2,'0')+rjust(str(event_start.hour),2,'0')+rjust(str(event_start.minute),2,'0')
            plt.savefig(out_path+event_string+'.png')
            plt.close()
                
        
        
            
    
            