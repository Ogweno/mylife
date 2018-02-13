from numpy import genfromtxt,zeros,c_,mean,arange
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace
from glob import glob


fgauges=glob('/Users/dmelgar/Tohoku2011/Tsunami/gauges/TI.*')
fout='/Users/dmelgar/Maule2010/wave_gauges/ancu.tsun'
time_epi=UTCDateTime('2011-03-11T05:46:23')
tcut=timedelta(seconds=7200)


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



for k in range(len(fgauges)):
    #Determine the kind of file
    temp=genfromtxt(fgauges[k],usecols=1,dtype='S')
    if temp[0]=='Mar':
        prs=genfromtxt(fgauges[k],usecols=4)
    else:
        prs=genfromtxt(fgauges[k],usecols=6)
    t=arange(0,len(prs)*15,15)
#plt.plot(t,prs-mean(prs),t,rad-mean(rad));plt.legend(['prs','rad']) ;plt.show()
#pick one
    data=prs.copy()
    #Interpolate to regular itnerval#Processing
    data=data-mean(data)
    fcorner=0.0001
    dt=t[2]-t[1]
    data_fil=highpass(data,fcorner,1./dt,2)
    plt.figure()
    plt.plot(t,data_fil)
    #Put in sac file
    station=fgauges[k].split('/')[-1].split('.')[1]
    plt.title(station)
    fout='/Users/dmelgar/Tohoku2011/Tsunami/gauges/'+station+'.tsun'
    st=Stream(Trace())
    st[0].stats.station=station
    st[0].stats.delta=dt
    st[0].data=data_fil
    st[0].data=st[0].data-st[0].data[0]
    st.write(fout,format='SAC')
    plt.show()