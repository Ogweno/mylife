from numpy import genfromtxt,mean
from glob import glob
from obspy import Stream,Trace,read
from obspy.core import UTCDateTime
from numpy.random import rand
from datetime import timedelta
from string import rjust

path=u'/Users/dmelgar/code/BBP/bbp/bbp_data/finished/rawdata/fake_nocal/_adjusted/'
sta=glob(path+'*acc.bbp')
noise_length=200
target_rate=40

#get noise
hypo_offset=timedelta(seconds=noise_length/2.)
noise=read(u'/Users/dmelgar/Data/BKS/fdsnws-dataselect_2016-07-13T01-35-52.mseed')
noise[0].decimate(factor=2)
noise_length=noise_length*40
#APPLY GAIN
noise[0].data=noise[0].data*100
noise[0].data=noise[0].data/214098.0
#noise[0].data=noise[0].data-mean(noise[0].data)
noise_start=int(rand(1)*0.8*len(noise[0].data))
n=noise.copy()
n[0].stats.delta=1./target_rate
st=Stream()
for k in range(len(sta)):
    #Pick a random pathc of noise
    noise_start=int(rand(1)*0.8*len(noise[0].data))
    n[0].data=noise[0].data[noise_start:noise_start+noise_length]
    #Bring noise_length/2 sample to zero to overlap with data
    noise_offset=n[0].data[noise_length/2]
    n[0].data=n[0].data-noise_offset
    #read station
    a=genfromtxt(sta[k])
    az=a[:,3]
    dt=a[1,0]-a[0,0]
    print dt
    tr=Trace()
    tr.data=az.copy()

    tr_len=len(tr.data)
    trout=Trace()
    trout.data=n[0].data.copy()
    #Add noise and data
    trout.data[noise_length/2:noise_length/2+tr_len]=trout.data[noise_length/2:noise_length/2+tr_len]+tr.data
    #Make integer
    trout.data=trout.data*1e6
    trout.data=trout.data.astype('i4')
    #trout.stats.mseed.encoding = u'INT32'
    t0=UTCDateTime('2016-09-07T14:42:26')
    trout.stats.starttime=t0-hypo_offset
    trout.stats.delta=dt
    trout.stats.station='ST'+rjust(str(k),2,'0')
    trout.stats.network='DM'
    trout.stats.channel='BNZ'
    print trout.id
    st+=trout
    

st.write('/Users/dmelgar/code/BBP/bbp/bbp_data/finished/rawdata/fake_nocal/_adjusted/fake_nocal_pfix.mseed',format='MSEED',encoding=11)    



