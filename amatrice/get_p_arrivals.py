from obspy.core import UTCDateTime
from obspy import read
from numpy import genfromtxt,where,arange,ones,zeros,array,tile,argmin
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

#time_epi=UTCDateTime('2016-08-24T01:36:32')
time_epi=UTCDateTime('2016-08-24T01:36:31.5')
epicenter=[13.214403,42.719490,8.0]
#epicenter=[13.24,42.70,8.1]

stations=['LSS','RM33','SPD','TERO','ASP','PTI','MNF','FEMA','FOS','TRE','SPM','NRC','AMT','GSA','PZI1']
Nsta=len(stations)
relative_ptimes=[5.831,65.7839,5.4050,67.8967,5.3254,6.5378,5.7208,63.1667,5.8395,5.9275,5.8311,5.825,5.66964,5.22,5.39]
rootpath='/Users/dmelgar/Amatrice2016/strong_motion/sac/'

lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[1,2])
station_catalogue=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta',usecols=[0],dtype='S')

#velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
velmod = TauPyModel(model='aci')

predictedP=9999*ones(len(stations))
#Get predicted arrivals
for ksta in range(len(stations)):
    st=read(rootpath+stations[ksta]+'.HNZ.sac')
    
    #Find coordinates
    i=where(station_catalogue==stations[ksta])[0]
    lon_sta=lonlat[i,0]
    lat_sta=lonlat[i,1]
    
                    
    deg=locations2degrees(lon_sta,lat_sta,epicenter[0],epicenter[1])
    arrivals = velmod.get_travel_times(source_depth_in_km=epicenter[2],distance_in_degree=deg,phase_list=['P','Pn','p'])
    
    #Observed P
    obsP=st[0].stats.starttime+relative_ptimes[ksta]
    
    #Determine P and S arrivals
    for kphase in range(len(arrivals)):
        if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
            if arrivals[kphase].time<predictedP[ksta]:
                predictedP[ksta]=arrivals[kphase].time
                
    print '%s: predicted=%s , observed=%s , difference=%.2f' % (stations[ksta],time_epi+predictedP[ksta],obsP,time_epi+predictedP[ksta]-obsP)
