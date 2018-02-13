from glob import glob
from obspy import read

#files=glob('/Users/dmelgar/Amatrice2016/strong_motion/sac/*HNE*')
#out='/Users/dmelgar/Amatrice2016/strong_motion/stations/latest.sta'
#f=open(out,'w')
#for k in range (len(files)):
#    st=read(files[k])
#    sta=files[0].split('/')[-1].split('.')[0]
#    lon=st[0].stats.longitude
#    lat=st[0].stats.latitude
#    line='%s\t%.6f\t%.6f\n' % (sta,lon,lat)
#    f.write(line)
#f.close()

#M5.3
files1=glob('/Users/dmelgar/Amatrice2016/M5.3/strong_motion/sac/*HNE*')
files2=glob('/Users/dmelgar/Amatrice2016/M5.3/strong_motion/sac/*HGE*')
out='/Users/dmelgar/Amatrice2016/M5.3/strong_motion/stations/latest.sta'
f=open(out,'w')
for k in range (len(files)):
    st=read(files[k])
    sta=files[0].split('/')[-1].split('.')[0]
    lon=st[0].stats.longitude
    lat=st[0].stats.latitude
    line='%s\t%.6f\t%.6f\n' % (sta,lon,lat)
    f.write(line)
f.close()