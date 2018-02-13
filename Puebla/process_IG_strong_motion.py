from obspy import read
from numpy import genfromtxt,where
from glob import glob
from obspy import UTCDateTime
from datetime import timedelta
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees



#read stations
path=u'/Users/dmelgar/Puebla2017/strong_motion/IG/'
allfiles=glob(path+'*.HLE.sac')

path_out=u'/Users/dmelgar/Puebla2017/strong_motion/sac/'
path_out_plot=u'/Users/dmelgar/Puebla2017/strong_motion/plots/'
origin=UTCDateTime('2017-09-19T18:14:38')
t1=timedelta(seconds=-20)
t2=timedelta(seconds=240)

#event location
hypo=[-98.643,18.313,45]

model = TauPyModel(model="iasp91")


#Read stations file
stafile=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/IG_stations.txt',usecols=1,dtype='S')
lons=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/IG_stations.txt',usecols=2)
lats=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/IG_stations.txt',usecols=3)

for ksta in range(len(allfiles)):

    sta=allfiles[ksta].split('.')[0].split('/')[-1]
    i=where(stafile==sta)[0]
    stalat=lats[i]
    stalon=lons[i]
    
    D=locations2degrees(hypo[1],hypo[0],stalat,stalon)
    arrivals = model.get_travel_times(source_depth_in_km=hypo[2],distance_in_degree=D,
                                  phase_list=["P", "p"])
    
    t0=timedelta(seconds=arrivals[0].time)        

    #let's do east
    component='HLE'
    files=glob(path+'*'+sta+'*'+component+'*')
    print sta
    for ktrace in range(len(files)):
        if ktrace==0:
            st=read(files[ktrace])
        else:
            st+=read(files[ktrace])
    st.merge(fill_value=0)
    st[0].trim(starttime=origin+t1+t0,endtime=origin+t2+t0)
    st[0].data=st[0].data
    e=st.copy()
    st.write(path_out+sta+'.'+component+'.sac',format='SAC')
        
        
    #let's do east
    component='HLN'
    files=glob(path+'*'+sta+'*'+component+'*')
    for ktrace in range(len(files)):
        if ktrace==0:
            st=read(files[ktrace])
        else:
            st+=read(files[ktrace])
    st.merge(fill_value=0)
    st[0].trim(starttime=origin+t1+t0,endtime=origin+t2+t0)
    st[0].data=st[0].data
    n=st.copy()
    st.write(path_out+sta+'.'+component+'.sac',format='SAC')
        
        
    #let's do east
    component='HLZ'
    files=glob(path+'*'+sta+'*'+component+'*')
    for ktrace in range(len(files)):
        if ktrace==0:
            st=read(files[ktrace])
        else:
            st+=read(files[ktrace])
    st.merge(fill_value=0)
    st[0].trim(starttime=origin+t1+t0,endtime=origin+t2+t0)
    st[0].data=st[0].data
    z=st.copy()
    st.write(path_out+sta+'.'+component+'.sac',format='SAC')
        

    #plot
        
    max_ampl=max([max(abs(e[0].data)),max(abs(n[0].data)),max(abs(z[0].data))])
    
    plt.figure(figsize=(12,6))
    
    plt.subplot(311)
    plt.plot(e[0].times(),e[0].data)
    plt.grid()
    plt.ylabel('East (m/s/s)')
    plt.ylim([-max_ampl,max_ampl])
    
    plt.subplot(312)
    plt.plot(n[0].times(),n[0].data)
    plt.grid()
    plt.ylabel('North (m/s/s)')
    plt.ylim([-max_ampl,max_ampl])
    
    plt.subplot(313)
    plt.plot(z[0].times(),z[0].data)
    plt.grid()
    plt.ylabel('Vert. (m/s/s)')
    plt.xlabel('Seconds')
    plt.ylim([-max_ampl,max_ampl])
    
    plt.suptitle(sta)
    
    plt.savefig(path_out_plot+sta+'.png')
    plt.close()
            
        