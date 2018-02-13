from numpy import genfromtxt,zeros,c_,where
from obspy.core import UTCDateTime
from obspy import Stream,Trace,read
from scipy.interpolate import interp1d


stations=genfromtxt('/Users/dmelgar/tidegauge_noise/station_info/stations.txt',usecols=1,dtype='S')
path_out='/Users/dmelgar/tidegauge_noise/data/'
months=['01','02','03','04','05','06','07','08','09','10','11','12']

# whatchu wanna do?
convert_months=False
paste_all=True
despike=True

#variable='tide'
variable='pressure'
#variable='temperature'



def highpass(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt



def process_tide_gauge(fin,fout,station):


    time=c_[genfromtxt(fin,usecols=0,dtype='S'),genfromtxt(fin,usecols=1,dtype='S')]
    t1=UTCDateTime(time[0,0]+time[0,1])
    #Make every time relative ins econds to t1
    t=zeros(len(time))
    for k in range(len(time)):
        t[k]=UTCDateTime(time[k,0]+time[k,1])-t1
    
    #read data
    data=genfromtxt(fin,usecols=2)
    
    #metadata?
    dt=int(t[1]-t[0])
    
    #Put in sac file
    st=Stream(Trace())
    st[0].stats.station=station
    st[0].stats.delta=dt
    st[0].stats.starttime=t1
    st[0].data=data
    st.write(fout,format='SAC')


if convert_months:
    for sta in stations:
        print sta
        if sta=='cres' or sta=='porf' or sta=='nspi':
            for kmonth in range(len(months)):
                if variable=='tide':
                    suffix=''
                elif variable=='pressure':
                    suffix='pres'
                elif variable=='temperature':
                    suffix='temp'
                try:
                    fin=path_out+sta+'/csv/'+sta+'_'+suffix+'_'+months[kmonth]+'_2017.csv'
                    fout=path_out+sta+'/sac/'+sta+'_'+suffix+'_'+months[kmonth]+'_2017.sac'
                    process_tide_gauge(fin,fout,sta)
                except:
                    print 'Error at station '+sta+' on month '+months[kmonth]
            
            
            
            
if paste_all:
    for sta in stations:
        if sta=='cres' or sta=='porf' or sta=='nspi':
            print sta
            for kmonth in range(len(months)):  
                if variable=='tide':
                    suffix=''
                elif variable=='pressure':
                    suffix='pres'
                elif variable=='temperature':
                    suffix='temp'
                
                if kmonth==0:
                    st=read(path_out+sta+'/sac/'+sta+'_'+suffix+'_'+months[kmonth]+'_2017.sac')
                elif kmonth!=1:
                    st+=read(path_out+sta+'/sac/'+sta+'_'+suffix+'_'+months[kmonth]+'_2017.sac')
    
            #st.merge()
            st.write(path_out+sta+'/sac/'+sta+'_'+suffix+'_2017.mseed',format='mseed')
            
            
            
if despike:
    for sta in stations:
        if sta=='cres' or sta=='porf' or sta=='nspi':
            print sta
            for kmonth in range(len(months)):  
                if variable=='tide':
                    suffix=''
                elif variable=='pressure':
                    suffix='pres'
                elif variable=='temperature':
                    suffix='temp'
                
                st=read(path_out+sta+'/sac/'+sta+'_'+suffix+'_'+'2017.mseed')

                for k in range(len(st)):
                    #Where are the NOT gaps
                    i=where(st[k].data!=1)[0]
                    #interpolate tot he whole thing
                    interpolant=interp1d(st[k].times()[i],st[k].data[i],bounds_error=False)
                    tinterp=st[k].times()
                    data_interp=interpolant(tinterp)
                    st[k].data=data_interp
    
            #st.merge()
            st.write(path_out+sta+'/sac/'+sta+'_'+suffix+'_2017.mseed',format='mseed')
