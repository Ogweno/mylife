from numpy import mean
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from datetime import timedelta
from obspy import read
from glob import glob


pathin='/Users/dmelgar/Iquique2014/tsunami/gauge_data/'
pathout='/Users/dmelgar/Iquique2014/tsunami/gauge_data/45/'
time_epi=UTCDateTime('2014-04-01T23:46:47')
tcut=timedelta(seconds=45*60)

files=glob(pathin+'*.sac')
for k in range(len(files)):
    sta=files[k].split('/')[-1].split('.')[0]
    st=read(files[k])
    st.trim(starttime=time_epi,endtime=time_epi+tcut)
    if sta=='aric':
        st[0].data=st[0].data-mean(st[0].data[0:10])
    plt.figure()
    plt.plot(st[0].times()/60,st[0].data)
    plt.grid()
    plt.xlabel('Minutes after OT')
    plt.ylabel('Sea surface height (m)')
    plt.title(st[0].stats.station)
    st.write(pathout+sta+'.tsun',format='SAC')
    plt.savefig(pathout+sta+'.png')
    plt.close()
