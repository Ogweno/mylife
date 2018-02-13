from glob import glob
from numpy import genfromtxt,where
from obspy import Stream,Trace,read,UTCDateTime
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees
from datetime import timedelta

folders_list=glob('/Users/dmelgar/Maule2010/strong_motion/cesmd/raw/*')
sta_names=genfromtxt('/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=0,dtype='S')
dist=genfromtxt('/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=3)
delta=genfromtxt('/Users/dmelgar/Maule2010/strong_motion/sm.sta',usecols=5)
out_folder='/Users/dmelgar/Maule2010/strong_motion/cesmd/proc/'
hypo=[-72.733002,-35.909000,35]
hypo_time=UTCDateTime('2010-02-2T06:34:11Z')


makesac=False
consolidate=True

def cosmos2sac(cosmos_file,sta_name=None):
    '''
    convert from cosmos format to sac
    '''
    
    from obspy import Stream,Trace,UTCDateTime
    from numpy import array,nan,isnan
    
    #Read file
    f=open(cosmos_file,'r')
    data=[]
    while True:
        line=f.readline()
        if line=='':
            break
        #Found the date/time of the first sample
        if 'DataSeries.FirstSampleTime.DateTime_txt' in line:
            date=line.split('=')[1].replace('"','').replace(';','')
            date=UTCDateTime(date)
        #Found sample rate
        if 'DataSeries.SamplesPerSecond_dbl' in line:
            dt=1./float(line.split()[2])
        #Found dt
        if 'RawSeries.SampleInterval_dbl' in line:
            dt=float(line.split()[2])
        #Azimuth and inclination
        if 'Sensor.Azimuth.Value_dbl' in line:
            if 'NULL' in line:
                azimuth=nan
            else:
                azimuth=float(line.split()[2])
        if 'Sensor.Inclination.Value_dbl' in line:
            if 'NULL' in line:
                inclination=nan
            else:
                inclination=float(line.split()[2])
        #Data gain
        if 'DataSeries.OrdinateUnits(1)_txt' in line:
            if 'g_standard/10' in line:
                gain=0.1*9.81
            elif 'cm/s/s' in line:
                gain=0.01
            else:
                print 'ERROR: unkown gain units '+line
        #Found start of data
        if '{' in line:
            while True:
                line=f.readline()
                if '}' in line: #end of data
                    break
                #Append to data object
                data.append(float(line))
    f.close()
    
    #Fix azimuth and inclination using file name
    if isnan(azimuth)==True or isnan(inclination)==True:
        if 'ch1' in cosmos_file.lower():
            azimuth=90
            inclination=90
        if 'ch2' in cosmos_file.lower():
            azimuth=0
            inclination=90
        if 'ch3' in cosmos_file.lower():
            azimuth=0
            inclination=0
             
    print azimuth,inclination                               
       
    #So what component was that??
    if azimuth==360:
        azimuth=0
    if azimuth>0 and azimuth<45:
        azimuth=0
    if azimuth>60 and azimuth<135:
        azimuth=90
    
    if azimuth==90:
        channel='HNE'
    if azimuth==0 and inclination==90:
        channel='HNN'
    if azimuth==0 and inclination==0:
        channel='HNZ'
                   
    data=array(data)*gain
    st=Stream(Trace())
    st[0].data=data
    st[0].stats.delta=dt
    st[0].stats.starttime=date
    st[0].stats.channel=channel
    st[0].stats.station=sta_name
    
    return st
    
    
if makesac:

    for k in range(len(folders_list)):
        
        print folders_list[k]
        
        files=glob(folders_list[k]+'/*cosm')
        station=folders_list[k].split('/')[-1]
        
        for kf in range(len(files)):
                
            st=cosmos2sac(files[kf],sta_name=station)
            
            fout=str(st[0].stats.starttime)+'.'+st[0].stats.station+'.'+st[0].stats.channel
            fout=folders_list[k]+'/'+fout+'.sac'
            
            st.write(fout,format='SAC')
            
            
            
if consolidate:
    
    velmod=TauPyModel(model='/Users/dmelgar/Maule2010/strong_motion/maule')
    
    for k in range(len(folders_list)):
        
        print folders_list[k]
        
        # East
        files=glob(folders_list[k]+'/*HNE.sac')
        e=Stream(Trace())
        for kf in range(len(files)):
            e+=read(files[kf])
        e.merge(fill_value='interpolate')
        
        # north
        files=glob(folders_list[k]+'/*HNN.sac')
        n=Stream(Trace())
        for kf in range(len(files)):
            n+=read(files[kf])
        n.merge(fill_value='interpolate')
        
        # East
        files=glob(folders_list[k]+'/*HNZ.sac')
        z=Stream(Trace())
        for kf in range(len(files)):
            z+=read(files[kf])
        z.merge(fill_value='interpolate')
        
        #get hypo to station distance
        i=where(sta_names==e[0].stats.station)[0]
        dist_in_km=dist[i]
        dist_in_degs=kilometer2degrees(dist_in_km)
        Ppaths=velmod.get_ray_paths(hypo[2],dist_in_degs,phase_list=['P','p'])
        first_sample_time=hypo_time+timedelta(seconds=Ppaths[0].path['time'][-1])+timedelta(seconds=delta[k])
        
        e[0].stats.starttime=first_sample_time
        n[0].stats.starttime=first_sample_time
        z[0].stats.starttime=first_sample_time
        
        eout=out_folder+e[0].stats.station+'.'+e[0].stats.channel+'.sac'
        nout=out_folder+n[0].stats.station+'.'+n[0].stats.channel+'.sac'
        zout=out_folder+z[0].stats.station+'.'+z[0].stats.channel+'.sac'
        
        e.write(eout,format='sac')
        n.write(nout,format='sac')
        z.write(zout,format='sac')
            