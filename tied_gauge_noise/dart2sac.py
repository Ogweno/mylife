from numpy import genfromtxt,array,where,arange
from obspy import Stream,Trace,UTCDateTime
from scipy.interpolate import interp1d

sta='46404'
mo='06'
fin='/Users/dmelgar/tidegauge_noise/DART/txt/'+sta+'_2017_'+mo+'.txt'
fout='/Users/dmelgar/tidegauge_noise/DART/sac/'+sta+'_2017_'+mo+'.sac'

dt1=15
dt2=15*60

dart=genfromtxt(fin)
eta=dart[:,-1]

#start time=
t1=UTCDateTime(str(dart[-1,0].astype('int'))+'-'+str(dart[-1,1].astype('int'))+'-'+str(dart[-1,2].astype('int'))+'T'+str(dart[-1,3].astype('int'))+':'+str(dart[-1,4].astype('int'))+':'+str(dart[-1,5].astype('int')))

#convert to UTC Date times
t=[]
t_string=[]
for k in range(len(dart)):
    tnow=UTCDateTime(str(dart[k,0].astype('int'))+'-'+str(dart[k,1].astype('int'))+'-'+str(dart[k,2].astype('int'))+'T'+str(dart[k,3].astype('int'))+':'+str(dart[k,4].astype('int'))+':'+str(dart[k,5].astype('int')))
    t_string.append(tnow)
    t.append(tnow-t1)
    
t=array(t)
t=t[::-1]
eta=eta[::-1]

#remove bad data
i=where(eta<9999)[0]
t=t[i]
eta=eta[i]
t0=t.min()
t1=t.max()
t_string=array(t_string)
t_string=t_string[::-1]
t_string=t_string[i]

tout=arange(t0,t1,dt1)
f=interp1d(t,eta)
etaout=f(tout)


#lowpass filter and decimate from 15s to 15mins
fcorner=1./dt2
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

etaout=lowpass(etaout,fcorner,1./dt1,2)
tout2=arange(t0,t1,dt2)
f=interp1d(tout,etaout)
etaout2=f(tout2)


#sac
st=Stream(Trace())
st[0].data=etaout2
st[0].stats.delta=dt2
st[0].stats.starttime=t_string[0]
st.write(fout,format='SAC')