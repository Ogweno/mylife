from numpy import genfromtxt,where,zeros,nan,ones
from glob import glob
from obspy.core.util.geodetics import gps2DistAzimuth
from matplotlib import pyplot as plt
from obspy import read
from obspy.core import UTCDateTime
from datetime import timedelta

lonepi=-122.3174
latepi=38.2118
time_epi=UTCDateTime('2014-08-24T10:20:44')
tplot=timedelta(seconds=100)
mul=1.5

pgd=genfromtxt('/Users/dmelgar/Napa2014/PGD/napa_test_nolatency.txt')
path='/Users/dmelgar/Napa2014/GPS/sac/'
lonlat=genfromtxt(u'/Users/dmelgar/Napa2014/unr_coords.txt',usecols=[1,2])
lon=lonlat[:,0]
lat=lonlat[:,1]
stas=genfromtxt(u'/Users/dmelgar/Napa2014/unr_coords.txt',usecols=0,dtype='S')
#Get lsit of files
filesn=glob(path+'*LXN.sac')
filese=glob(path+'*LXE.sac')
#Initalize
d=zeros(len(filese)) #epicentral distances
#Loop and plot
dmin=[]
dmax=0
plt.figure()
f,axarr=plt.subplots(1,2)
axe=axarr[1]
axn=axarr[0]
for k in range(len(filese)):
    current_sta=filese[k].split("/")[-1].split(".")[0].upper()
    i=where(current_sta==stas)[0]
    try:
        d,az,baz=gps2DistAzimuth(latepi,lonepi,lat[i],lon[i])
        d=d/1000
        dmin=min([dmin,d])
        dmax=max([dmax,d])
    except:
        d=nan
    #Read data
    stn=read(filesn[k])
    ste=read(filese[k])
    #Trim
    stn.trim(starttime=time_epi,endtime=time_epi+tplot,pad=True,fill_value=0)
    ste.trim(starttime=time_epi,endtime=time_epi+tplot,pad=True,fill_value=0)
    #Self Normalize
    stn[0].data=stn[0].data/max([stn[0].data.max(),-stn[0].data.min()])
    ste[0].data=ste[0].data/max([ste[0].data.max(),-ste[0].data.min()])
    dplot=ones(ste[0].times().shape)*d
    #Plot
    axn.plot(stn[0].times(),stn[0].data*mul+dplot,'k')
    axe.plot(ste[0].times(),ste[0].data*mul+dplot,'k')
axn.set_title('North')
axe.set_title('East')
axn.set_ylim(dmin-5,75)
axe.set_ylim(dmin-5,75)
axn.grid()
axe.grid()
axn.set_xlabel('Seconds after OT')
axe.set_xlabel('Seconds after OT')
axn.set_ylabel('Epicentral distance (km)')
axe.yaxis.set_ticklabels([])
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.05, hspace=0)

fig, ax1 = plt.subplots()
ax1.scatter(pgd[:,1],pgd[:,2])
ax1.set_xlabel('Seconds after OT')
ax1.set_xlim(0,100)
ax1.set_ylabel('Mw', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.scatter(pgd[:,1], pgd[:,3],marker='+', c='r')
ax2.set_ylabel('No. stations', color='r')
ax2.set_ylim(0,50)
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.set_xlim(0,100)
plt.show()
    

