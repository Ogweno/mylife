from glob import glob
from numpy import genfromtxt,where,array,float64,diff,where,mean
from obspy import Stream,Trace,read
from os import remove
from obspy.core import UTCDateTime
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees
from datetime import timedelta
from mudpy.forward import lowpass
from matplotlib import pyplot as plt
from shutil import move

path='/Users/dmelgar/Napa2014/GPS/sac/'
gps_files=glob(path+'GPS/difkin_????')
time_epi=UTCDateTime('2014-08-24T10:20:44')
t0=UTCDateTime('2010-04-04T00:00:00')
epicenter=array([-115.287,32.259,10]) 


make_pgd=True


        
if make_pgd:
    pathout='/Users/dmelgar/PGD/GPS/sac/Napa2014/'
    tcut=timedelta(minutes=10)
    files=glob(path+'*LXN.sac')
    for k in range(len(files)):
        sta=files[k].split('/')[-1].split('.')[0]
        n=read(path+sta+'.LXN.sac')
        e=read(path+sta+'.LXE.sac')
        u=read(path+sta+'.LXZ.sac')
        #Trim
        n[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        e[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        u[0].trim(starttime=time_epi,endtime=time_epi+tcut)
        #Remove mean of first ten seconds
        n[0].data=n[0].data-mean(n[0].data[0:16])
        e[0].data=e[0].data-mean(e[0].data[0:16])
        u[0].data=u[0].data-mean(u[0].data[0:16])
        #Write to file
        n.write(pathout+sta+'.LXN.sac',format='SAC')
        e.write(pathout+sta+'.LXE.sac',format='SAC')
        u.write(pathout+sta+'.LXZ.sac',format='SAC')