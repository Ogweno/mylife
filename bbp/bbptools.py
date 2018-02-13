
def read_srf(srf_file):

    from numpy import array,r_,zeros,argmin,where
    '''
    Parse SRF file (oh joy)
    '''

    f=open(srf_file)
    
    begin_read=False
    line1=f.readline()
    kpoint=0
    stf_all=[]
    while True:
        if begin_read==False:
            line1=f.readline()
        if 'POINTS' in line1:
            points=int(line1.split()[1])
            begin_read=True
            xyz=zeros((points,3))
            slip=zeros(points)
            tinit=zeros(points)
            rise_time=zeros(points)
        if begin_read:
            #Read infor about this point source location, type, etc
            line1=f.readline()
            if line1 == '':
                break
            #Assign latlon points
            xyz[kpoint,0]=float(line1.split()[0])
            xyz[kpoint,1]=float(line1.split()[1])
            xyz[kpoint,2]=float(line1.split()[2])
            tinit[kpoint]=float(line1.split()[6])
            dt=float(line1.split()[7])
            line2=f.readline()
            slip[kpoint]=float(line2.split()[1])
            kpoint+=1
            #Number of STF points
            Nstf=int(line2.split()[2])
            #Format only allows 6 pts epr line, so how many lines do we need to read?
            if Nstf % 6==0:
                Nread=Nstf/6
            else:
                Nread=Nstf/6+1
            #Read stf for this point source
            stf=array([])
            for kread in range(Nread):
                line=f.readline()
                stf=r_[stf,array(line.split()).astype('float')]
            stf_all.append(stf)
            rise_time[kpoint-1]=len(stf)*dt
    #Find hypocenter
    i=where(tinit>0)[0]
    imin=argmin(tinit[i])
    hypocenter=xyz[i[imin],:]
    f.close()
    
    return xyz,slip,tinit,stf_all,rise_time,hypocenter
  
  
def get_simulation_intensity(sim_path):
    '''
    Get PGA, PGV and Rjb for a simulation
    '''  
  
    from glob import glob
    from numpy import genfromtxt,zeros,linspace,ones,r_
    from scipy.integrate import cumtrapz
    from pyproj import Geod   

    stations=genfromtxt(glob(sim_path+'/param_files/*.stl')[0],usecols=2,dtype='S')
    lonlat=genfromtxt(glob(sim_path+'/param_files/*.stl')[0],usecols=[0,1])
    
    #WHere's the src file
    src_file=glob(sim_path+'/param_files/*.src')[0]
    
    #Parse src file for magnitude, center coords and fault dims
    f=open(src_file)
    while True:
        line=f.readline()
        if 'MAGNITUDE' in line:
            M=float(line.split('=')[1])
        if 'FAULT_LENGTH' in line:
            length=float(line.split('=')[1])
        if 'FAULT_WIDTH' in line:
            width=float(line.split('=')[1])
        if 'STRIKE' in line:
            strike=float(line.split('=')[1])
        if 'LAT_TOP' in line:
            lat_top=float(line.split('=')[1])
        if 'LON_TOP' in line:
            lon_top=float(line.split('=')[1])
            f.close()
            break

    #Make fault trace for calculating Rjb
    g=Geod(ellps='WGS84')
    dist=linspace(0.1,length/2,200)*1000
    lon1,lat1,foo=g.fwd(lon_top*ones(len(dist)),lat_top*ones(len(dist)),strike*ones(len(dist)),dist)
    lon2,lat2,foo=g.fwd(lon_top*ones(len(dist)),lat_top*ones(len(dist)),strike*ones(len(dist))-180,dist)
    lon_trace=r_[lon1,lon2]
    lat_trace=r_[lat1,lat2]
    
    #Initalize
    pga=zeros(len(stations))
    pgv=zeros(len(stations))
    pgd=zeros(len(stations))
    Rjb=zeros(len(stations))
    
    #Get ID#
    files=glob(sim_path+'/*.bbp')
    ID=files[0].split('/')[-1].split('.')[0]
    
    for kstation in range(len(stations)):
        
        current_station=stations[kstation]
        
        v=genfromtxt(sim_path+'/'+ID+'.'+current_station+'.vel.bbp')
        a=genfromtxt(sim_path+'/'+ID+'.'+current_station+'.acc.bbp')
        
        #Get PGA in g and PGV in cm/s
        pgv[kstation]=v[:,1:3].max()
        pga[kstation]=a[:,1:3].max()/981
        
        #Integrate displacememnts
        dx=cumtrapz(v[:,2],v[:,0],initial=0)
        dy=cumtrapz(v[:,1],v[:,0],initial=0)
        dz=cumtrapz(v[:,3],v[:,0],initial=0)
        d=(dx**2+dy**2+dz**2)**0.5
        
        # Get pgd in cm
        pgd[kstation]=d.max() 
        
        #Get Rjb
        sta_lon=lonlat[kstation,0]*ones(len(lon_trace))
        sta_lat=lonlat[kstation,1]*ones(len(lon_trace))
        foo,bar,trace_dist=g.inv(sta_lon,sta_lat,lon_trace,lat_trace)
        Rjb[kstation]=trace_dist.min()/1000
        
    return pga,pgv,pgd,Rjb,M

        
def adjust_to_ptime(path_to_sim,stl_file,srf_file):  
    '''
    Fix pre-event jitter in LF seismograms
    '''
 
    from numpy import genfromtxt,where,diff,sign,zeros,r_,c_,savetxt
    from glob import glob
    import os
    import shutil
    
    #Get stationd ata
    stations=genfromtxt(stl_file,usecols=2,dtype='S')
    lonlat=genfromtxt(stl_file,usecols=[0,1])
    
    #Loop over stations
    for sta in stations:
        
        #LF file names
        lf_north_file=glob(path_to_sim+'*'+sta+'*lf-resamp.000')[0]
        lf_east_file=glob(path_to_sim+'*'+sta+'*lf-resamp.090')[0]
        lf_up_file=glob(path_to_sim+'*'+sta+'*lf-resamp.ver')[0]
        
        #HF file names
        hf_north_file=glob(path_to_sim+'*'+sta+'*hf-resamp.000')[0]
        hf_east_file=glob(path_to_sim+'*'+sta+'*hf-resamp.090')[0]
        hf_up_file=glob(path_to_sim+'*'+sta+'*hf-resamp.ver')[0]   
        
        #Read LF files
        tlf,lf_north=read_bbp_seismogram(lf_north_file)
        tlf,lf_east=read_bbp_seismogram(lf_east_file)
        tlf,lf_up=read_bbp_seismogram(lf_up_file)
        
        #Read HF files
        thf,hf_north=read_bbp_seismogram(hf_north_file)
        thf,hf_east=read_bbp_seismogram(hf_east_file)
        thf,hf_up=read_bbp_seismogram(hf_up_file)
        
        #Get p time for this station
        i=where(stations==sta)[0]
        xyz,slip,tinit,stf_all,rise_time,hypocenter=read_srf(srf_file)
        ptime,stime=arrivals(hypocenter,lonlat[i,0],lonlat[i,1])
        
        #Data sampling interval
        dt=tlf[1]-tlf[0]

        
        #Find zero crossing in LF data tot he LEFT of the ptime
        i=where(tlf<ptime)
        zero_crossings = where(diff(sign(lf_up[i])))[0]
        
        #Keep closest one to ptime
        zero_crossing=zero_crossings[-1]
        tshift=ptime-tlf[zero_crossing]
        
        #What is thift in dt units?
        tshift/dt
        tshift=int(round(tshift/dt))
        
        #Add this many zeros and chop off htis many samples from the end
        z=zeros(tshift)
        lf_north=r_[z,lf_north[0:-tshift]]
        lf_east=r_[z,lf_east[0:-tshift]]
        lf_up=r_[z,lf_up[0:-tshift]]
        
        #Mute everything before this time
        lf_north[i]=0
        lf_east[i]=0
        lf_up[i]=0
        lf_north[i[-1]+1]=0
        lf_east[i[-1]+1]=0
        lf_up[i[-1]+1]=0
        
        #Add hf and lf
        npts=len(lf_north)
        N=lf_north+hf_north[0:npts]
        E=lf_east+hf_east[0:npts]
        Z=lf_up+hf_up[0:npts]
        T=tlf.copy()
        
        # Check if folder exists, make it, write new bbp files
        folder=path_to_sim+'_adjusted/'
        if not os.path.exists(folder):
            os.makedirs(folder)
            
        header='ptime adjusted BB sim\ntime(sec)      N-S(cm/s)      E-W(cm/s)      U-D(cm/s)'
        savetxt(folder+sta+'.acc.bbp',c_[T,N,E,Z],fmt='%.8e\t%.8e\t%.8e\t%.8e',header=header)
        
        
                  
                                
        
def read_bbp_seismogram(f,return_dt=False):
    '''
    Read intermediate output files from BBP runs
    '''
    
    from numpy import zeros,arange
    
    seis=open(f,'r')
    #Skip first line
    seis.readline()
    #Read header data
    header=seis.readline()
    #How many points in this seismogram
    npts=int(header.split()[0])
    #Sampling interval in seconds
    dt=float(header.split()[1])
    
    #Read file and build output array and time vector
    t=arange(0,npts*dt,dt)
    out=zeros(npts)
    k=0
    while True:
        
        line=seis.readline()
        if line=='':
            break
        
        #How many points in this current data line?
        points_per_line=len(line.split())
        #Read line and convert to float
        line_floats=zeros(points_per_line)
        for j in range(points_per_line):
            line_floats[j]=float(line.split()[j])
        out[k:k+points_per_line]=line_floats
        
        #Update copunter and let's go again
        k+=points_per_line
    if return_dt:
        return t,out,dt
    else:
        return t,out
    
    
          

def hypo_position(src_in,hypo_fraction=0.05):
    '''
    determine the along strike position of the hypocenter
    '''
    
    f=open(src_in)
    while True:
        line=f.readline()
        if 'FAULT_LENGTH' in line:
            length=float(line.split('=')[1])
            f.close()
            break
    hypo_pos=hypo_fraction*length-length/2
    
    return hypo_pos


def arrivals(hypocenter,station_lon,station_lat):
    '''
    Get P and S arrival times
    '''
    
    from obspy.taup import TauPyModel
    from pyproj import Geod
    from numpy import rad2deg
    
    #Station to hypo distance
    g=Geod(ellps='WGS84')
    az,baz,dist=g.inv(hypocenter[0],hypocenter[1],station_lon,station_lat)
    #Convert distance from m to degrees and km
    dist_deg=rad2deg(dist/6371e3)
    print dist/1000
    #Calculate theoretical arrivals
    mod=TauPyModel(model='Nocal')
    arrivals=mod.get_travel_times(source_depth_in_km=hypocenter[2], distance_in_degree=dist_deg, phase_list=('P','p','S','s'))
    print arrivals
    
    #Parse arrivals
    ptime=1e6
    stime=1e6
    for k in range(len(arrivals)):
        if arrivals[k].name=='p' or arrivals[k].name=='P':
            ptime=min(arrivals[k].time,ptime)
        if arrivals[k].name=='s' or arrivals[k].name=='S':
            stime=min(arrivals[k].time,stime)
    print ptime
                   
    return ptime,stime
    
    

def paste_seismograms(lf_file,hf_file,filter_corner=1.0,filter_order=2,two_pass=False,seismogram_type='v',dt_synth=0.01):   
    '''
    filter and paste hf and lf seismograms
    '''
    
    from scipy.interpolate import interp1d
    from scipy.integrate import cumtrapz
    from numpy import r_,diff,arange
    
    #read data
    tlf,lf,dt_lf=read_bbp_seismogram(lf_file,return_dt=True)
    thf,hf,dt_hf=read_bbp_seismogram(hf_file,return_dt=True)
    
    #filter data
    lf_filter=filter_waveform(lf,filter_corner,1./dt_lf,filter_order,two_pass=two_pass,filter_type='low')
    hf_filter=filter_waveform(hf,filter_corner,1./dt_hf,filter_order,two_pass=two_pass,filter_type='high')

    # time interpolation vector
    t=arange(0,min(tlf[-1],thf[-1]),dt_synth)

    #sum and diff/integrate depending on type
    if seismogram_type=='a': #make acceleration
        flf=interp1d(tlf, lf_filter, kind='cubic', bounds_error=False)
        fhf=interp1d(thf, hf_filter, kind='linear', bounds_error=False)
        lf_interp=flf(t)
        hf_interp=fhf(t)
        seis=lf_interp+hf_interp    
        seis=r_[0,diff(seis)/dt_synth]       
                  
    elif seismogram_type=='v': #keep as velocity
        flf=interp1d(tlf, lf_filter, kind='linear', bounds_error=False)
        fhf=interp1d(thf, hf_filter, kind='linear', bounds_error=False)
        lf_interp=flf(t)
        hf_interp=fhf(t)
        seis=lf_interp+hf_interp 
    elif seismogram_type=='d': #integrate to displacement
        flf=interp1d(tlf, lf_filter, kind='linear', bounds_error=False)
        fhf=interp1d(thf, hf_filter, kind='linear', bounds_error=False)
        lf_interp=flf(t)
        hf_interp=fhf(t)
        seis=lf_interp+hf_interp  
        seis=cumtrapz(seis, t, dt_synth, initial=0) 
    
         
    return t,seis
    
    



def filter_waveform(data,fcorner,fsample,order,two_pass=True,filter_type='low'):
    '''
    Make a highpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt,lfilter
    from numpy import size,array
    
    fnyquist=fsample/2
    if filter_type=='low':
        b, a = butter(order, array(fcorner)/(fnyquist),'lowpass')
    elif filter_type=='high':
        b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
        
    if two_pass==True:
        data_filt=filtfilt(b,a,data)
    else:
        data_filt=lfilter(b,a,data)
    
    return data_filt
    