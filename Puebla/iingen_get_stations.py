from glob import glob
from numpy import zeros

files=glob(u'/Users/dmelgar/Puebla2017/strong_motion/II/Datos_2017-09-19-1814-40M7.1Axochiapan,Mor/*')
fout=open(u'/Users/dmelgar/Puebla2017/strong_motion/station_info/II_stations.txt','w')
lat=zeros(len(files))
lon=zeros(len(files))

for k in range(len(files)):
    if k==0:
        sta=[files[k].split('/')[-1][0:4]]
    else:
        sta.append(files[k].split('/')[-1][0:4])
    with open(files[k]) as fp:
        for i, line in enumerate(fp):
            if i == 22:
                lat[k]=float(line.split()[5])
            elif i == 23:
                lon[k]=float(line.split()[1])
            elif i > 23:
                break
                
for k in range(len(files)):
    line='%s\t%.5f\t%.5f\n' % (sta[k],-lon[k],lat[k])
    print line
    fout.write(line)
fout.close()