from obspy import read
from glob import glob

files=glob(u'/Users/dmelgar/Puebla2017/strong_motion/raw/*e170919.s13')
path_out=u'/Users/dmelgar/Puebla2017/strong_motion/sac_noproc/'



for k in range(len(files)):
    st=read(files[k])
    
    print st
    
    sta=st[0].stats.station
    
    chan=st[0].stats.channel
    st[0].write(path_out+sta+'.'+chan+'.sac',format='SAC')
    
    chan=st[1].stats.channel
    st[1].write(path_out+sta+'.'+chan+'.sac',format='SAC')
    
    chan=st[2].stats.channel
    st[2].write(path_out+sta+'.'+chan+'.sac',format='SAC')
    
    