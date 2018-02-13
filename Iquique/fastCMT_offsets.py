from numpy import genfromtxt,array,r_,ones,mean,arange,sqrt,savetxt,c_
from obspy import read
from run_filt import RunningMedian as med
from matplotlib import pyplot as plt

stations=genfromtxt('/Users/dmelgar/PGD/station_info/iquique.sta',usecols=0,dtype='S')
lon=genfromtxt('/Users/dmelgar/PGD/station_info/iquique.sta',usecols=1)
lat=genfromtxt('/Users/dmelgar/PGD/station_info/iquique.sta',usecols=2)
path= u'/Users/dmelgar/PGD/GPS/sac/Iquique2014/'
fout="/Users/dmelgar/Iquique2014/fastCMT/offsets.txt"
pathfig="/Users/dmelgar/Iquique2014/fastCMT/plots/"
tmed=30
tavg=10
twait=180
thresh=0.0
#Intialize
nout=[]
eout=[]
uout=[]
lonout=[]
latout=[]
keep=[]
for k in range(len(stations)):
    try:
        e=read(path+stations[k]+'.LXE.sac')
        n=read(path+stations[k]+'.LXN.sac')
        u=read(path+stations[k]+'.LXZ.sac')
        #Calcualte moving average
        em=e.copy()
        nm=n.copy()
        um=u.copy()
        em[0].data=r_[ones(tmed-1)*e[0].data[0],array(med(e[0].data,tmed))]
        nm[0].data=r_[ones(tmed-1)*n[0].data[0],array(med(n[0].data,tmed))]
        um[0].data=r_[ones(tmed-1)*u[0].data[0],array(med(u[0].data,tmed))]
        eoffset=mean(em[0].data[twait:twait+tavg])
        noffset=mean(nm[0].data[twait:twait+tavg])
        uoffset=mean(um[0].data[twait:twait+tavg])
        toffset=arange(twait,twait+tavg)
        #Make plot
        #plt.figure(figsize=(15,8))
        #plt.subplot(311)
        #plt.plot(e[0].times(),e[0].data,em[0].times(),em[0].data)
        #plt.ylabel('East(m)')
        #plt.plot(toffset,ones(len(toffset))*eoffset,c='r',lw=2)
        #plt.xlim([0,300])
        #plt.grid()
        #plt.subplot(312)
        #plt.plot(n[0].times(),n[0].data,nm[0].times(),nm[0].data)
        #plt.ylabel('North(m)')
        #plt.plot(toffset,ones(len(toffset))*noffset,c='r',lw=2)
        #plt.xlim([0,300])
        #plt.grid()
        #plt.subplot(313)
        #plt.plot(u[0].times(),u[0].data,um[0].times(),um[0].data)
        #plt.ylabel('Up(m)')
        #plt.plot(toffset,ones(len(toffset))*uoffset,c='r',lw=2)
        #plt.xlim([0,300])
        #plt.grid()
        #plt.savefig(pathfig+stations[k]+'.median.pdf')
        #plt.close()
        # Now save things
        print sqrt(noffset**2+eoffset**2)
        if sqrt(noffset**2+eoffset**2)>thresh:
            eout.append(eoffset)
            nout.append(noffset)
            uout.append(uoffset)
            lonout.append(lon[k])
            latout.append(lat[k])
    except:
        print 'Station '+stations[k]+' not found'
#save to file
keep=array(keep)
lonout=array(lonout)
latout=array(latout)
eout=array(eout)
nout=array(nout)
uout=array(uout)
f=open(fout,'w')
for k in range(len(lon)):
    line='%s\t%12.6f\t%12.6f\t%8.4f\t%8.4f\t%8.4f\n' %(stations[k],lon[k],lat[k],nout[k],eout[k],uout[k])
    f.write(line)