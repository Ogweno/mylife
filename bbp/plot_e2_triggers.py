from matplotlib import pyplot as plt
from numpy import genfromtxt,array,where
from obspy.core import UTCDateTime
from pyproj import Geod
from numpy import rad2deg
from obspy.taup import TauPyModel
from obspy import read
import bbptools

#Parse triggers
f=open('/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/triggers_last.txt')
stations=[]
triggers=[]
while True:
    line=f.readline()
    if line=='':
        break
    if 'BNZ' in line:
        stations.append(line.split()[5])
        triggers.append(line.split()[11])


#Station coordinatesinsa
all_stations=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/param_files/fake_nocal.stl',usecols=[2],dtype='S')
lonlat=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/param_files/fake_nocal.stl',usecols=[0,1])

#Station miniSEED
st=read(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/rawdata/fake_nocal/_adjusted/fake_nocal_pfix.mseed')

#Station to hypo distance
g=Geod(ellps='WGS84')
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/large_eew_m5.0_frac0.5.srf')

xl=[-0.02,13]
t0=UTCDateTime('2016-09-07T14:42:26')
fig, axarr = plt.subplots(20, 1,figsize=(8,16))
for k in range(len(stations)):   
    a=genfromtxt(u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/fake_nocal/111.'+stations[k].lower()+'.vel.bbp')
    trigger_time=UTCDateTime(triggers[k])-t0
    
    
    #Find current station
    i=where(all_stations==stations[k].lower())[0]

    az,baz,dist=g.inv(hypocenter[0],hypocenter[1],lonlat[i,0],lonlat[i,1])
    #Convert distance from m to degrees and km
    dist_deg=rad2deg(dist/6371e3)
    dist_km=dist/1000
    
    #Calculate theoretical arrivals
    mod=TauPyModel(model='Nocal')
    arrivals=mod.get_travel_times(source_depth_in_km=hypocenter[2], distance_in_degree=dist_deg, phase_list=('P','p','S','s'))
    
    #Parse arrivals
    for ka in range(len(arrivals)):
        if arrivals[ka].name=='p' or arrivals[ka].name=='P':
            ptime=arrivals[ka].time
        if arrivals[ka].name=='s' or arrivals[ka].name=='S':
            stime=arrivals[ka].time
    
    #Find trace for plotting
    kname=0
    while True:
        if st[kname].stats.station.upper()==stations[k].upper():
            tr=st[kname].copy()
            tr.data=tr.data/1e7
            break
        kname+=1
    
    ax=axarr[k]
    
    ax.plot(tr.times()-100,tr.data)
    yl=array([tr.data.min(),tr.data.max()])/10
    ax.plot([trigger_time,trigger_time],[-1e10,1e10],'r',lw=2)
    ax.plot([ptime,ptime],[-1e10,1e10],'k',lw=2)
    ax.set_ylim(yl)
    ax.set_xlim(xl)
    ax.get_yaxis().set_visible(False)
    if k<19:
        ax.xaxis.set_ticklabels([])
        
    if k==19:
        ax.set_xlabel('Seconds after OT')
    
plt.show()
    
        