from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,where,arange,ones,zeros,array,tile,argmin
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

time_epi=UTCDateTime('2016-08-24T01:36:32')
time_offsets=[-0.6,-0.4,-0.2,0,0.2,0.4,0.6]

stations=['LSS','RM33','SPD','TERO','ASP','PTI','MNF','FEMA','FOS','TRE','SPM','NRC','AMT']
Nsta=len(stations)
relative_ptimes=[5.831,65.7839,5.4050,67.8967,5.3254,6.5378,5.7208,63.1667,5.8395,5.9275,5.8311,5.825,5.66964]
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/sac/'

lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[1,2])
station_catalogue=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[0],dtype='S')

#velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
velmod = TauPyModel()

lon_grid=arange(13.17,13.37,0.01)
lat_grid=arange(42.68,42.75,0.01)
z_grid=arange(2,10,0.5)
Npts=len(lon_grid)*len(lat_grid)*len(z_grid)*len(time_offsets)*len(stations)

predictedP=ones((Nsta,Npts))*9999
arrival_times=zeros((Nsta,Npts))
error=zeros((Nsta,Npts))

lon_out=zeros(Npts)
lat_out=zeros(Npts)
z_out=zeros(Npts)
offsets_out=zeros(Npts)



#Get absolute times

for ksta in range(len(stations)):
    st=read(rootpath+stations[ksta]+'.HNZ.sac')
    
    #Find coordinates
    i=where(station_catalogue==stations[ksta])[0]
    lon_sta=lonlat[i,0]
    lat_sta=lonlat[i,1]
    
    #Loop over sources, get travel times
    ksource=0
    for ktime in range(len(time_offsets)):
        print 'Working on station %d and time offset %d' %(ksta,ktime)
        for klon in range(len(lon_grid)):
            for klat in range(len(lat_grid)):
                for kz in range(len(z_grid)):
                    
                    
                    #Apply time offsets
                    arrival_times[ksta,ksource]=st[0].stats.starttime+relative_ptimes[ksta]-time_epi+time_offsets[ktime]
                    
                    lon_out[ksource]=lon_grid[klon]
                    lat_out[ksource]=lat_grid[klat]
                    z_out[ksource]=z_grid[kz]
                    offsets_out[ksource]=time_offsets[ktime]
                    
                    deg=locations2degrees(lon_sta,lat_sta,lon_grid[klon],lat_grid[klat])
                    arrivals = velmod.get_travel_times(source_depth_in_km=z_grid[kz],distance_in_degree=deg,phase_list=['P','Pn','p'])
                    
                    #Determine P and S arrivals
                    for kphase in range(len(arrivals)):
                        if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
                            if arrivals[kphase].time<predictedP[ksta,ksource]:
                                predictedP[ksta,ksource]=arrivals[kphase].time
                    error[ksta,ksource]=predictedP[ksta,ksource]-arrival_times[ksta,ksource]
                                
                    ksource+=1
                
                
#calculate L2 error for each source point
weights=array([1,1,1,1,1,1,1,1,1,1,1,1,1])
weights=tile(weights,(Npts,1)).T
arrival_times=array(arrival_times)
arrival_times=tile(arrival_times,(Npts,1)).T
L2=(((predictedP-arrival_times)**2).sum(axis=0))**0.5


for k in range(len(z_grid)):
    j=where(z_out==z_grid[k])[0]
    i=argmin(L2[j])
    min_point=j[i]
    print 'lon=%.4f\tlat=%.4f\tz=%.4f\tt_offset=%.2f\tL2=%.6f' % (lon_out[min_point] , lat_out[min_point] , z_out[min_point],offsets_out[min_point],L2[min_point])
