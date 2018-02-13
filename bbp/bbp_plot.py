from matplotlib.colors import LinearSegmentedColormap as lsc

cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.20, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.20, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.075, 1, 1),
                  (0.20, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
whitejet = lsc('whitejet',cdict,256)

def plot_1sta(sim_path,staname,hypocenter,stafile,tlims=[0,30]):
    '''
    Plot one station and it's predicted arrivals from ray tracing
    '''
    from matplotlib import pyplot as plt
    from obspy.taup import TauPyModel
    from numpy import genfromtxt,where,rad2deg
    from glob import glob
    from pyproj import Geod

    #load data
    afile=glob(sim_path+'*.'+staname+'.acc*')[0]
    a=genfromtxt(afile)
    vfile=glob(sim_path+'*.'+staname+'.vel*')[0]
    v=genfromtxt(vfile)

    #Station coordinates
    all_stations=genfromtxt(sim_path+'param_files/'+stafile,usecols=[2],dtype='S')
    i=where(all_stations==staname)[0]
    lonlat=genfromtxt(sim_path+'param_files/'+stafile,usecols=[0,1])[i][0]
    
    #Station to hypo distance
    g=Geod(ellps='WGS84')
    az,baz,dist=g.inv(hypocenter[0],hypocenter[1],lonlat[0],lonlat[1])
    #Convert distance from m to degrees
    dist_deg=rad2deg(dist/6371e3)
    
    #Calculate theoretical arrivals
    mod=TauPyModel(model='Nocal')
    arrivals=mod.get_travel_times(source_depth_in_km=hypocenter[2], distance_in_degree=dist_deg, phase_list=('P','p','S','s'))
    
    #Parse arrivals
    for k in range(len(arrivals)):
        if arrivals[k].name=='p' or arrivals[k].name=='P':
            ptime=arrivals[k].time
        if arrivals[k].name=='s' or arrivals[k].name=='S':
            stime=arrivals[k].time
            
    #Pplot the thingamajig
    
    plt.figure(figsize=(12,6))
    
    plt.subplot(321)
    plt.plot(a[:,0],a[:,1],'k')
    plt.ylabel('North')
    plt.title('Accel. (cm/s/s)')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)
    
    plt.subplot(323)
    plt.plot(a[:,0],a[:,2],'k')
    plt.ylabel('East')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)
    
    plt.subplot(325)
    plt.plot(a[:,0],a[:,3],'k')
    plt.ylabel('Up')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)
            
            
    plt.subplot(322)
    plt.plot(v[:,0],v[:,1],'k')
    plt.ylabel('North')
    plt.title('Vel. (cm/s)')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)
    
    plt.subplot(324)
    plt.plot(v[:,0],v[:,2],'k')
    plt.ylabel('East')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)
    
    plt.subplot(326)
    plt.plot(v[:,0],v[:,3],'k')
    plt.ylabel('Up')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)
    plt.xlim(tlims)

    plt.show()
    
    
def plot_pwave(sim_path,staname,stafile,tlims=[0,30],ylims=[-2,2],filter=True,corner=1./13):
    '''
    Plot one station and it's predicted arrivals from ray tracing
    '''
    from matplotlib import pyplot as plt
    from obspy.taup import TauPyModel
    from numpy import genfromtxt,where,rad2deg,ones,log10,argmin
    from scipy.integrate import cumtrapz
    from glob import glob
    from pyproj import Geod
    import bbptools
    from mudpy.forward import highpass


    #Get hypocenter
    srf_file=glob(sim_path+'*.srf')[0]
    xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)
    i=argmin(tinit)
    hypocenter=xyz[i,:]

    #load data
    vfile=glob(sim_path+'*.'+staname+'.vel*')[0]
    afile=glob(sim_path+'*.'+staname+'.acc*')[0]
    a=genfromtxt(afile)
    v=genfromtxt(vfile)
    t=v[:,0]
    vz=v[:,3]
    dt=t[1]-t[0]
    accz=a[:,3]
    ta=a[:,0]
    
    #Get displacement
    dz=cumtrapz(vz,t,initial=0)
    
    #Filter?
    if filter==True:
        dfil=highpass(dz,corner,1/dt,4)
        dz=dfil.copy()

    #Station coordinates
    all_stations=genfromtxt(sim_path+'param_files/'+stafile,usecols=[2],dtype='S')
    i=where(all_stations==staname)[0]
    lonlat=genfromtxt(sim_path+'param_files/'+stafile,usecols=[0,1])[i][0]
    
    #Station to hypo distance
    g=Geod(ellps='WGS84')
    az,baz,dist=g.inv(hypocenter[0],hypocenter[1],lonlat[0],lonlat[1])
    #Convert distance from m to degrees and km
    dist_deg=rad2deg(dist/6371e3)
    dist_km=dist/1000
    
    #Calculate theoretical arrivals
    mod=TauPyModel(model='Nocal')
    arrivals=mod.get_travel_times(source_depth_in_km=hypocenter[2], distance_in_degree=dist_deg, phase_list=('P','p','S','s'))
    
    #Parse arrivals
    for k in range(len(arrivals)):
        if arrivals[k].name=='p' or arrivals[k].name=='P':
            ptime=arrivals[k].time
        if arrivals[k].name=='s' or arrivals[k].name=='S':
            stime=arrivals[k].time
  
    #Get expected value of Pd
    metadata=glob(sim_path+'param_files/*.src')[0]          
    f=open(metadata,'r')
    line=f.readline()
    mag=float(line.split()[-1])
    f.close()
    Pd=10**((mag-5.39-1.38*log10(dist_km))/1.23)
    
    #Assign ylims
    ylims=[-2*Pd,2*Pd]
                                                                                                            
    #Plot the thingamajig  
    plt.figure(figsize=(9,4.5))
    plt.subplot(212)
    
    plt.plot(t,dz,'k')
    plt.plot(t,Pd*ones(len(t)),'--')
    plt.plot(t,-Pd*ones(len(t)),'--')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)    
    plt.ylabel('Vertical (cm)')
    plt.xlabel('Seconds after OT')
    
    plt.xlim(tlims)
    plt.ylim(ylims)
    
    plt.subplot(211)
        
    plt.plot(ta,accz,'k')
    plt.scatter(ptime,0,marker='|',c='b',s=100,lw=2)
    plt.scatter(stime,0,marker='|',c='r',s=100,lw=2)    
    plt.ylabel('Vertical (cm/s/s)')
    
    plt.xlim(tlims)
    
    plt.subplots_adjust(left=0.12,right=0.97,bottom=0.12,top=0.96,hspace=0.12)
    plt.show()
    


def plot_allP(sim_path,sim_list):
    '''
    Extract all P-values and scatter plot
    '''
    from glob import glob
    from numpy import genfromtxt,where,log10,r_,ones,array,logspace,argmin
    import bbptools
    from scipy.integrate import cumtrapz
    from pyproj import Geod
    from matplotlib import pyplot as plt
    
    
    #Initalize
    observed_Pd=[]
    predicted_Pd=[]
    M=array([])
    R=array([])
    
    #Loop over sim_list
    for ksim in sim_list:
        print ksim
        stations=genfromtxt(glob(sim_path+ksim+'/param_files/*.stl')[0],usecols=2,dtype='S')
        lonlat=genfromtxt(glob(sim_path+ksim+'/param_files/*.stl')[0],usecols=[0,1])
        
        #Get ID#
        files=glob(sim_path+ksim+'/*.bbp')
        ID=files[0].split('/')[-1].split('.')[0]
        
        for kstation in range(len(stations)):
            
            current_station=stations[kstation]
            
            v=genfromtxt(sim_path+ksim+'/'+ID+'.'+current_station+'.vel.bbp')
            t=v[:,0]
            vz=v[:,3]
            #Get displacement
            dz=cumtrapz(vz,t,initial=0)
            
            #Read meta data get magnitude and hypocenter
            if kstation==0:
                #Get hypocenter
                srf_file=glob(sim_path+ksim+'/*.srf')[0]
                xyz,slip,tinit,stf_all=bbptools.read_srf(srf_file)
                i=argmin(tinit)
                hypocenter=xyz[i,:]

                #get magnitude
                metadata=glob(sim_path+ksim+'/param_files/*.src')[0]          
                f=open(metadata,'r')
                line=f.readline()
                mag=float(line.split()[-1])
                f.close()
                M=r_[M,mag*ones(len(stations))]
            
            #Get P,S arrivals
            ptime,stime=arrivals(hypocenter,lonlat[kstation,0],lonlat[kstation,1])
            
            #Get  observed Pd
            ptime=ptime-1
            i=where((t>=ptime) & (t<=ptime+3))[0]
            Pd=max(abs(dz[i]))
            observed_Pd.append(Pd)
            
            #Get predicted Pd
            
            #Station to hypo distance
            g=Geod(ellps='WGS84')
            az,baz,dist=g.inv(hypocenter[0],hypocenter[1],lonlat[kstation,0],lonlat[kstation,1])
            #Convert distance from m to degrees and km
            dist_km=dist/1000
            R=r_[R,dist_km]
                
            
            #Finally, actually get Pd
            Pd=10**((mag-5.39-1.38*log10(dist_km))/1.23)      
            predicted_Pd.append(Pd)  
            
    observed_Pd=array(observed_Pd)
    predicted_Pd=array(predicted_Pd)
            
    #Make plot
    plt.figure()    
    i=where(M==4.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DC143C')
    i=where(M==5.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#32CD32')
    i=where(M==5.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#0000CD')
    i=where(M==6.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DAA520')
    i=where(M==6.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#9932CC')
    i=where(M==7.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#202020')

    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    #Reference lines
    Rref=logspace(0,2)
    Pd45=10**((4.5-5.39-1.38*log10(Rref))/1.23)
    Pd50=10**((5.0-5.39-1.38*log10(Rref))/1.23) 
    Pd55=10**((5.5-5.39-1.38*log10(Rref))/1.23) 
    Pd60=10**((6.0-5.39-1.38*log10(Rref))/1.23) 
    Pd65=10**((6.5-5.39-1.38*log10(Rref))/1.23) 
    Pd70=10**((7.0-5.39-1.38*log10(Rref))/1.23) 
    
    plt.plot(Rref,Pd45,'--',lw=2,c='#DC143C')
    plt.plot(Rref,Pd50,'--',lw=2,c='#32CD32')
    plt.plot(Rref,Pd55,'--',lw=2,c='#0000CD')
    plt.plot(Rref,Pd60,'--',lw=2,c='#DAA520')
    plt.plot(Rref,Pd65,'--',lw=2,c='#9932CC')
    plt.plot(Rref,Pd70,'--',lw=2,c='#202020')
    
    ax.set_xlim([1e-4,1e1])
    ax.set_xlim([1e0,1e2])
    plt.xlabel('Distance (km)')
    plt.ylabel('Pd (cm)')
    
    plt.legend(['M4.5','M5.0','M5.5','M6.0','M6.5','7.0'],loc=3)
    
    
    
    plt.figure(figsize=(10,6))
    
    plt.subplot(231)
    i=where(M==4.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DC143C')
    plt.plot(Rref,Pd45,'--',lw=2,c='#DC143C')
    plt.legend(['M4.5'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([1e0,1e2])
    plt.ylabel('Pd (cm)')
    
    plt.subplot(232)
    i=where(M==5.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#32CD32')
    plt.plot(Rref,Pd50,'--',lw=2,c='#32CD32')
    plt.legend(['M5.0'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([1e0,1e2])
    
    plt.subplot(233)
    i=where(M==5.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#0000CD')
    plt.plot(Rref,Pd55,'--',lw=2,c='#0000CD')
    plt.legend(['M5.5'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([1e0,1e2])
    
    plt.subplot(234)
    i=where(M==6.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DAA520')
    plt.plot(Rref,Pd60,'--',lw=2,c='#DAA520')
    plt.legend(['M6.0'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlabel('Distance (km)')
    plt.ylabel('Pd (cm)')
    ax.set_xlim([1e0,1e2])
    
    plt.subplot(235)
    i=where(M==6.5)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#9932CC')
    plt.plot(Rref,Pd65,'--',lw=2,c='#9932CC')
    plt.legend(['M6.5'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlabel('Distance (km)')
    ax.set_xlim([1e0,1e2])
    
    plt.subplot(236)
    i=where(M==7.0)[0]
    plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#202020')
    plt.plot(Rref,Pd70,'--',lw=2,c='#202020')
    plt.legend(['M7.0'],frameon=False,loc=3)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlabel('Distance (km)')
    ax.set_xlim([1e0,1e2])
    
    plt.show()
    
    
            
def plot_srf(srf_file,src_file,contour_interval=1.0,highlight=3.0,UTM_zone='10S',figsize=(8,4),xlim=[-20,20],save=False,fout=None):
    '''
    plot srf file
    '''
    
    from matplotlib import pyplot as plt
    import bbptools
    from pyproj import Proj,Geod
    from numpy import linspace,ones,r_,zeros,argmin,where,arange
    from mudpy.viewFQ import pqlx
    import matplotlib.tri as tri
    from matplotlib import ticker
    
    #Parse src file
    f=open(src_file)
    while True:
        line=f.readline()
        if 'FAULT_LENGTH' in line:
            length=float(line.split('=')[1])*1000
        if 'STRIKE' in line:
            strike=float(line.split('=')[1])
            break
    
    xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)
    
    #Project fault coordiantes to strike direction
    p=Proj(ellps='WGS84',proj='utm',zone=UTM_zone)
    g=Geod(ellps='WGS84')
    
    #First make straight line to be projected onto
    npts=10000
    lon_orig=hypocenter[0]*ones(npts)
    lat_orig=hypocenter[1]*ones(npts)
    strike=r_[(strike-180)*ones(npts/2),strike*ones(npts/2)]
    L=r_[linspace(-1.2*length,-0.1,npts/2),linspace(0.1,1.1*length,npts/2)]
    lon_line,lat_line,foo=g.fwd(lon_orig,lat_orig,strike,L)
    lat_line=linspace(xyz[:,1].min()-1.0,xyz[:,1].max()+1.0,10000)
    
    #Convert everything to UTM
    x_line,y_line=p(lon_line,lat_line)
    x_hypo,y_hypo=p(hypocenter[0],hypocenter[1])
    x_source,y_source=p(xyz[:,0],xyz[:,1])
    
    # Get closest coordinates on line
    x_proj=zeros(len(x_source))
    y_proj=zeros(len(x_source))
    for k in range(len(x_source)):
        dist=((x_source[k]-x_line)**2+(y_source[k]-y_line)**2)**0.5
        i=argmin(dist)
        x_proj[k]=x_line[i]
        y_proj[k]=y_line[i]
    z_proj=-xyz[:,2]  
    
    #hypocenter projection
    i=where(tinit>0)[0]
    imin=argmin(tinit[i])
    x_hypo_proj=x_proj[i[imin]]
    y_hypo_proj=y_proj[i[imin]]
    z_hypo_proj=-xyz[i[imin],2]
    
    #Convert to along strike and down dip distances
    strike_dist=((x_proj-x_hypo_proj)**2+(y_proj-y_hypo_proj)**2)**0.5
    i=where(y_proj>y_hypo_proj)[0]
    strike_dist[i]=-strike_dist[i]
    #strike_dist[i]=strike_dist[i]
    strike_dist/=1000
    
    dip_dist=(((z_proj-z_hypo_proj)*1000)**2)**0.5
    i=where(z_proj<z_hypo_proj)[0]
    dip_dist[i]=-dip_dist[i]
    dip_dist/=1000
    
    # Make plot
    triang = tri.Triangulation(strike_dist, dip_dist)
    xl=xlim
    
    plt.figure(figsize=figsize)
    plt.subplot(211)
    plt.tripcolor(triang, slip, shading='flat', cmap=whitejet)
    plt.ylabel('Down dip (km)')
    cb=plt.colorbar(shrink=0.85)
    cb.set_label('Slip (cm)')
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    
    yl=[dip_dist.min(),dip_dist.max()]
    plt.ylim(yl)
    plt.xlim(xl)
    
    #Rupture onset contours
    plt.tricontour(triang,tinit,levels=arange(0,tinit.max()+contour_interval,contour_interval),colors='k',linewidths=1)
    plt.tricontour(triang,tinit,levels=[highlight],colors='k',linewidths=3)
    
    plt.scatter(0,0,c='w',marker='*',lw=1,s=250)
    
    
    ax=plt.gca()
    #ax.set_aspect('equal', adjustable='box')
    
    
    
    
    plt.subplot(212)
    plt.tripcolor(triang, rise_time, shading='flat', cmap=plt.cm.magma_r)
    plt.ylabel('Down dip (km)')
    plt.xlabel('Along Strike (km)')
    cb=plt.colorbar(shrink=0.85)
    cb.set_label('Rise time (s)')
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    
    yl=[dip_dist.min(),dip_dist.max()]
    plt.ylim(yl)
    plt.xlim(xl)
    
    #Rupture onset contours
    plt.tricontour(triang,tinit,levels=arange(0,tinit.max()+contour_interval,contour_interval),colors='k',linewidths=1)
    plt.tricontour(triang,tinit,levels=[highlight],colors='k',linewidths=3)
    
    plt.scatter(0,0,c='w',marker='*',lw=1,s=250)
    
    ax=plt.gca()
    #ax.set_aspect('equal', adjustable='box')
    plt.xticks(rotation=0)
    plt.subplots_adjust(hspace=0,wspace=0,top=0.96,bottom=0.11,left=0.1,right=0.98)
    
    plt.show()
    
    if save==True:
        plt.savefig(fout)
    
    
              
       
    
    
    
    
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
    print dist
    dist_deg=rad2deg(dist/6371e3)
    
    #Calculate theoretical arrivals
    mod=TauPyModel(model='Nocal')
    arrivals=mod.get_travel_times(source_depth_in_km=hypocenter[2], distance_in_degree=dist_deg, phase_list=('P','p','S','s'))
    
    #Parse arrivals
    for k in range(len(arrivals)):
        if arrivals[k].name=='p' or arrivals[k].name=='P':
            ptime=arrivals[k].time
        if arrivals[k].name=='s' or arrivals[k].name=='S':
            stime=arrivals[k].time
            
    return ptime,stime