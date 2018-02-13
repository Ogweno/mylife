from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from mtspec import mtspec
from numpy import mean


client = Client("IRIS")
#stations=['M18B','M10B','G34B','M13B','FS20B','FS13B']
stations=['M10B','G34B','M13B','FS20B','FS13B']
path='/Users/dmelgar/tidegauge_noise/APG/'

t1 = UTCDateTime("2013-05-01T00:00:00.000")
t2 = UTCDateTime("2013-06-01T00:00:00.000")

for k in range(len(stations)):

    sta=stations[k]
    print sta
    
    st = client.get_waveforms("7D", sta, "", "HDH", t1, t2)
    
    #got too mnay traces?
    st.merge()
    
    #gain goes from counts to Pa
    gain=99.7254
    st[0].data=st[0].data/gain
    st[0].data=st[0].data/(9.81*1029) # Pa -> m 
    
    fcorner=1./(60) #1 minute corner
    
    def lowpass(data,fcorner,fsample,order):
        '''
        Make a lowpass zero phase filter
        '''
        from scipy.signal import butter,filtfilt
        from numpy import size,array
        
        fnyquist=fsample/2
        b, a = butter(order, array(fcorner)/(fnyquist),'lowpass')
        data_filt=filtfilt(b,a,data)
        return data_filt
        
    #Anti alias filter with one minute corner
    st[0].data=lowpass(st[0].data,fcorner,1./st[0].stats.delta,2)
    
    #decimate to 1s sampling
    st[0].decimate(factor=125,no_filter=True)
    
    #decimate to 1 minute sampling
    st[0].decimate(factor=60,no_filter=True)
    
    st.write(path+sta+'_2013_05.sac',format='sac')



##Psd
#Tw=5
#Ntapers=8
#psd, f, jackknife, _, _ = mtspec(data=st[0].data, delta=st[0].stats.delta, time_bandwidth=Tw,number_of_tapers=Ntapers, nfft=st[0].stats.npts, statistics=True)