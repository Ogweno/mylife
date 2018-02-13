from obspy import read
from numpy import genfromtxt
from glob import glob
from obspy import UTCDateTime
from datetime import timedelta
from matplotlib import pyplot as plt

#read stations
net=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/all_stations.txt',usecols=0,dtype='S')
sta=genfromtxt('/Users/dmelgar/Puebla2017/strong_motion/station_info/all_stations.txt',usecols=1,dtype='S')
path=u'/Users/dmelgar/Puebla2017/strong_motion/raw/'
path_out=u'/Users/dmelgar/Puebla2017/strong_motion/sac/'
path_out_plot=u'/Users/dmelgar/Puebla2017/strong_motion/plots/'
origin=UTCDateTime('2017-09-19T18:14:38')
t1=timedelta(seconds=-20)
t2=timedelta(seconds=240)

for ksta in range(len(sta)):
    if net[ksta]=='VM':
        #let's do east
        plot=0
        component='HLE'
        files=glob(path+'*'+sta[ksta]+'*'+component+'*')
        if len(files)>0: #there is data
            plot+=1
            print sta[ksta]
            for ktrace in range(len(files)):
                if ktrace==0:
                    st=read(files[ktrace])
                else:
                    st+=read(files[ktrace])
            st.merge(fill_value=0)
            st[0].trim(starttime=origin+t1,endtime=origin+t2)
            st[0].data=st[0].data/100.
            e=st.copy()
            st.write(path_out+sta[ksta]+'.'+component+'.sac',format='SAC')
            
            
        #let's do east
        component='HLN'
        files=glob(path+'*'+sta[ksta]+'*'+component+'*')
        if len(files)>0: #there is data
            print sta[ksta]
            for ktrace in range(len(files)):
                if ktrace==0:
                    st=read(files[ktrace])
                else:
                    st+=read(files[ktrace])
            st.merge(fill_value=0)
            st[0].trim(starttime=origin+t1,endtime=origin+t2)
            st[0].data=st[0].data/100.
            n=st.copy()
            st.write(path_out+sta[ksta]+'.'+component+'.sac',format='SAC')
            
            
        #let's do east
        component='HLZ'
        files=glob(path+'*'+sta[ksta]+'*'+component+'*')
        if len(files)>0: #there is data
            print sta[ksta]
            for ktrace in range(len(files)):
                if ktrace==0:
                    st=read(files[ktrace])
                else:
                    st+=read(files[ktrace])
            st.merge(fill_value=0)
            st[0].trim(starttime=origin+t1,endtime=origin+t2)
            st[0].data=st[0].data/100.
            z=st.copy()
            st.write(path_out+sta[ksta]+'.'+component+'.sac',format='SAC')
            
        if plot>0:
            
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
            
            plt.suptitle(sta[ksta])
            
            plt.savefig(path_out_plot+sta[ksta]+'.png')
            plt.close()
            
        