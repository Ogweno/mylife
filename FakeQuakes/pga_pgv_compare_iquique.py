from obspy import read,UTCDateTime
from numpy import mean,where,log,log10,r_,genfromtxt,ones,zeros
from mudpy.forward import lowpass as bandpass
from mtspec import mtspec
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.integrate import cumtrapz
from pyproj import Geod



#stations=genfromtxt('/Volumes/Illapel/FQ/iquique/data/station_info/sm_5sta.gflist',usecols=0,dtype='S')
#lonlat=genfromtxt('/Volumes/Illapel/FQ/iquique/data/station_info/sm_5sta.gflist',usecols=[1,2])

stations=genfromtxt('/Volumes/Illapel/FQ/iquique/data/station_info/sm.gflist',usecols=0,dtype='S')
lonlat=genfromtxt('/Volumes/Illapel/FQ/iquique/data/station_info/sm.gflist',usecols=[1,2])

time_epi=UTCDateTime('2014-04-01T23:46:47Z')
epicenter=[-70.82,-19.64,30]

sim_path=u'/Volumes/Illapel/FQ/iquique/output/waveforms/iquique.000000/'
data_path=u'/Users/dmelgar/Iquique2014/SAC/PROC/'
tcut=120
fcorner=[1./10,24]
vmin=-7 ; vmax=0
cmap_spec=plt.cm.jet

mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14

#get epicentral dsitances
p=Geod(ellps='WGS84')
az,baz,d=p.inv(lonlat[:,0],lonlat[:,1],ones(len(lonlat))*epicenter[0],ones(len(lonlat))*epicenter[1])
d=d/1000.

pga_obs=zeros(len(stations))
pgv_obs=zeros(len(stations))
pga_fq=zeros(len(stations))
pgv_fq=zeros(len(stations))

for k in range(len(stations)):
    
    sta=stations[k]
    print sta

    #Read fakequakes
    nfq=read(sim_path+sta+'.bb.HNN.sac')
    efq=read(sim_path+sta+'.bb.HNE.sac')
    zfq=read(sim_path+sta+'.bb.HNZ.sac')

    #cut
    nfq.trim(starttime=time_epi,endtime=time_epi+tcut)
    efq.trim(starttime=time_epi,endtime=time_epi+tcut)
    zfq.trim(starttime=time_epi,endtime=time_epi+tcut)

    #Read real data
    n=read(data_path+sta+'.HNN.sac')
    e=read(data_path+sta+'.HNE.sac')
    z=read(data_path+sta+'.HNZ.sac')
    
    #Cut real data, remove pre-event mean
    n.trim(starttime=time_epi,endtime=time_epi+tcut)
    e.trim(starttime=time_epi,endtime=time_epi+tcut)
    z.trim(starttime=time_epi,endtime=time_epi+tcut)
    n[0].data=n[0].data-mean(n[0].data[0:1000])
    e[0].data=e[0].data-mean(e[0].data[0:1000])
    z[0].data=z[0].data-mean(z[0].data[0:1000])
    
    #get pga
    pga_obs[k]=r_[n[0].data,z[0].data,e[0].data].max()
    pga_fq[k]=r_[nfq[0].data,zfq[0].data,efq[0].data].max()
    
    #band pass filter both to same frequency band and integrate
    fsample=1./n[0].stats.delta
    order=4
    n[0].data=bandpass(n[0].data,fcorner,fsample,order,zerophase=True)
    e[0].data=bandpass(e[0].data,fcorner,fsample,order,zerophase=True)
    z[0].data=bandpass(z[0].data,fcorner,fsample,order,zerophase=True)
    
    fsample=1./nfq[0].stats.delta
    nfq[0].data=bandpass(nfq[0].data,fcorner,fsample,order,zerophase=True)
    efq[0].data=bandpass(efq[0].data,fcorner,fsample,order,zerophase=True)
    zfq[0].data=bandpass(zfq[0].data,fcorner,fsample,order,zerophase=True)
    
    #integrate
    n[0].data=cumtrapz(n[0].data,n[0].times(),initial=0)
    e[0].data=cumtrapz(e[0].data,e[0].times(),initial=0)
    z[0].data=cumtrapz(z[0].data,z[0].times(),initial=0)
    nfq[0].data=cumtrapz(nfq[0].data,nfq[0].times(),initial=0)
    efq[0].data=cumtrapz(efq[0].data,efq[0].times(),initial=0)
    zfq[0].data=cumtrapz(zfq[0].data,zfq[0].times(),initial=0)
    
    #get pgv
    pgv_obs[k]=r_[n[0].data,z[0].data,e[0].data].max()
    pgv_fq[k]=r_[nfq[0].data,zfq[0].data,efq[0].data].max()    


    
    
lims=[0.1,10]
fig, axarr = plt.subplots(1, 2,figsize=(12,4)) 
ax=axarr[0]
ax.set_yscale('log')
ax.set_xscale('log')
ax.scatter(pga_obs,pga_fq)
ax.plot(lims,lims)
ax.set_xlabel('Observed PGA')
ax.set_ylabel('FQ PGA')
for k in range(len(stations)):
    ax.annotate(s=stations[k],xy=(pga_obs[k],pga_fq[k]),fontsize=10)
ax.set_xlim(lims)
ax.set_ylim(lims)

lims=[0.01,1]
ax=axarr[1]
ax.set_yscale('log')
ax.set_xscale('log')
ax.scatter(pgv_obs,pgv_fq)
ax.plot(lims,lims)
ax.set_xlabel('Observed PGV')
ax.set_ylabel('FQ PGV')
for k in range(len(stations)):
    ax.annotate(s=stations[k],xy=(pgv_obs[k],pgv_fq[k]),fontsize=10)
ax.set_xlim(lims)
ax.set_ylim(lims)


plt.subplots_adjust(bottom=0.2,wspace=0.4)
plt.show()
