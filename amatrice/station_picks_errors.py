from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,where,arange,ones,zeros,array,tile,argmin
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

time_epi=UTCDateTime('2016-08-24T01:36:32')


stations=['LSS','RM33','SPD','TERO','ASP','PTI','MNF','FEMA','FOS','TRE','SPM','NRC','AMT']
Nsta=len(stations)
relative_ptimes=[5.831,65.7839,5.4050,67.8967,5.3254,6.5378,5.7208,63.1667,5.8395,5.9275,5.8311,5.825,5.66964]
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/sac/'

lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[1,2])
station_catalogue=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[0],dtype='S')

#velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
velmod = TauPyModel()

time_offset=0.0
lon_epi=13.24
lat_epi=42.7#42.71
z_epi=8.0



#Get absolute times
for ksta in range(len(stations)):
    predictedP=9999
    
    st=read(rootpath+stations[ksta]+'.HNZ.sac')
    
    #Find coordinates
    i=where(station_catalogue==stations[ksta])[0]
    lon_sta=lonlat[i,0]
    lat_sta=lonlat[i,1]
    
    #Apply time offsets
    arrival_time=st[0].stats.starttime+relative_ptimes[ksta]-time_epi+time_offset
    
    deg=locations2degrees(lon_sta,lat_sta,lon_epi,lat_epi)
    arrivals = velmod.get_travel_times(source_depth_in_km=z_epi,distance_in_degree=deg,phase_list=['P','Pn','p'])
    
    #Determine P and S arrivals
    for kphase in range(len(arrivals)):
        if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
            if arrivals[kphase].time<predictedP:
                predictedP=arrivals[kphase].time
    error=predictedP-arrival_time
    
    print 'Station %s, observed=%.4f, predicted=%.4f, error=%.4f' %(stations[ksta],arrival_time,predictedP,error)
                                