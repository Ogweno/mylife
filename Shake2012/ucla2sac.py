#Read in UCLA files and conert to SAC

from obspy import Stream,Trace
from glob import glob
from numpy import genfromtxt
from datetime import datetime

path_out='/Users/dmelgar/Shaketable 2012/Accel/sac/'
meta_file='/Users/dmelgar/Shaketable 2012/Accel/meta_split.txt'
files=glob('/Users/dmelgar/Shaketable 2012/Accel/ascii/*UCLA2.txt')
#Read metadata
gain=genfromtxt(meta_file,usecols=16)*9.8
level=genfromtxt(meta_file,usecols=8,dtype='S')
quadrant=genfromtxt(meta_file,usecols=9,dtype='S')
channel=genfromtxt(meta_file,usecols=7,dtype='S')
#Loop over files (This is gonna be slow dear)
for k in range(len(files)):
    print files[k]
    f=genfromtxt(files[k])
    t0=datetime.utcfromtimestamp(f[0,0])
    #Loop over channels
    for j in range(0,81):
        #Define station name
        if level[j]=='-1': #It's on shake table platten
            sta='PLAT'
        else:
            sta='L'+level[j]+'Q'+quadrant[j]
        #Define channel name and polarity multiplier
        if channel[j]=='X':
            chan='HNE'
            polarity=1
        elif channel[j]=='-X':
            cahn='HNE'
            polarity=-1
        elif channel[j]=='Y':
            chan='HNN'
            polarity=1
        elif channel[j]=='-Y':
            chan='HNN'
            polarity=-1
        elif channel[j]=='Z':
            chan='HNZ'
            polarity=1
        elif channel[j]=='-Z':
            chan='HNZ'
            polarity=-1
        #Read data into stream
        st=Stream(Trace())
        st[0].data=(f[:,j+1]*polarity)*gain[j]
        st[0].stats.delta=0.005
        st[0].stats.starttime=t0
        st[0].stats.station=sta
        st[0].stats.channel=chan
        #Write to file
        st.write(path_out+t0.strftime("%Y%m%d%H%M%S")+'.'+sta+'.'+chan+'.sac',format='SAC')
        
        
    
    