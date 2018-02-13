from gmpe_tools import bssa14
from numpy import genfromtxt,exp,log,ones,logspace,zeros,r_,array,where
from matplotlib import pyplot as plt
from pyproj import Geod
from obspy import read
from resp_spec import responseSpectrum
from scipy.integrate import cumtrapz
from matplotlib import rcParams

rcParams.update({'font.size': 12})

rootpath='/Users/dmelgar/Amatrice2016/M6.6/strong_motion/sac/'

lonlat=genfromtxt('/Users/dmelgar/Amatrice2016/M6.6/strong_motion/stations/latest.sta',usecols=[1,2])
stations=genfromtxt('/Users/dmelgar/Amatrice2016/M6.6/strong_motion/stations/latest.sta',usecols=[0],dtype='S')

#Source
epicenter=[13.088,42.855]
lonepi=ones(len(lonlat))*epicenter[0]
latepi=ones(len(lonlat))*epicenter[1]

#Rjb poly
Rjb_poly=genfromtxt('/Users/dmelgar/Amatrice2016/strong_motion/Rjb.txt')
Rjb_poly[:,0]=epicenter[0]
Rjb_poly[:,1]=epicenter[1]



fc=1./20
def highpass_filter(data,fcorner,fsample,order):
    '''
    Make a lowpass zero phase filter
    '''
    from numpy import size,array
    from scipy.signal import butter,filtfilt
    fnyquist=fsample/2
    b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
    data_filt=filtfilt(b,a,data)
    return data_filt



#loop over stations
kwrite=0
for k in range(len(lonlat)):
        
        print k
        
        n=read(rootpath+stations[k]+'.HNN.sac')
        e=read(rootpath+stations[k]+'.HNE.sac')
        nmax=n[0].data.max()
        emax=e[0].data.max()
        
        
        #Get spectral ordinates
        Sa03_n=responseSpectrum(n[0].stats.delta,n[0].data/9.81,[1./0.3], oscDamping=0.05)
        Sa03_e=responseSpectrum(e[0].stats.delta,e[0].data/9.81,[1./0.3], oscDamping=0.05)
        Sa10_n=responseSpectrum(n[0].stats.delta,n[0].data/9.81,[1./1.0], oscDamping=0.05)
        Sa10_e=responseSpectrum(e[0].stats.delta,e[0].data/9.81,[1./1.0], oscDamping=0.05)
        Sa30_n=responseSpectrum(n[0].stats.delta,n[0].data/9.81,[1./3.0], oscDamping=0.05)
        Sa30_e=responseSpectrum(e[0].stats.delta,e[0].data/9.81,[1./3.0], oscDamping=0.05)
    
        #get PGV
        ve=cumtrapz(e[0].data,e[0].times())
        ve=highpass_filter(ve,fc,1./e[0].stats.delta,2)
        vemax=ve.max()
        vn=cumtrapz(n[0].data,n[0].times())
        vn=highpass_filter(vn,fc,1./n[0].stats.delta,2)
        vnmax=vn.max()
    
        #get distances
        g=Geod(ellps='WGS84')
        az,baz,dist=g.inv(Rjb_poly[:,0],Rjb_poly[:,1],ones(len(Rjb_poly))*lonlat[k,0],ones(len(Rjb_poly))*lonlat[k,1])
        dist=min(dist)/1000.
        
        #Put in arrays
        if k==0:
            pga=array([max(emax,nmax)])
            pgv=array([max(vemax,vnmax)])
            Sa03=array([max(Sa03_n[0],Sa03_e[0])])
            Sa10=array([max(Sa10_n[0],Sa10_e[0])])
            Sa30=array([max(Sa30_n[0],Sa30_e[0])])
            Rjb_obs=array([dist])
            if lonlat[k,1]>epicenter[1]:
                n_or_s=array([1])
            else:
                n_or_s=array([-1])
        else:
            pga=r_[pga,array([max(emax,nmax)])]
            pgv=r_[pgv,array([max(vemax,vnmax)])]
            Sa03=r_[Sa03,array([max(Sa03_e[0],Sa03_n[0])])]
            Sa10=r_[Sa10,array([max(Sa10_e[0],Sa10_n[0])])]
            Sa30=r_[Sa30,array([max(Sa30_e[0],Sa30_n[0])])]
            Rjb_obs=r_[Rjb_obs,array([dist])]
            if lonlat[k,1]>epicenter[1]:
                n_or_s=r_[n_or_s,array([1])]
            else:
                n_or_s=r_[n_or_s,array([-1])]

#Get italy GMPE
italy=genfromtxt('/Users/dmelgar/Amatrice2016/GMPE/ShakeMap.Table')
i=where(italy[:,1]==6.3)[0]
Ritaly=italy[i,0]
pga_italy=10**(italy[i,6])/9.81
pgv_italy=10**(italy[i,5])*100
Sa03_italy=10**(italy[i,4])/9.81
Sa10_italy=10**(italy[i,3])/9.81
Sa30_italy=10**(italy[i,2])/9.81

#Scale to stupid units
pga=pga/9.81
pgv=pgv*100


#Get theoretical GMPE
npts=100
M=6.6*ones(npts)
Rjb=logspace(0,2.5,npts)
U=zeros(npts)
SS=zeros(npts)
RS=zeros(npts)
NS=ones(npts)
vs30=760*ones(npts)

#get intensities
pga_bssa14,stdev_pga=bssa14(M, Rjb, vs30,U=U,SS=SS,RS=RS,NS=NS,Z1=None,intensity_measure='PGA',italy=True)
pgv_bssa14,stdev_pgv=bssa14(M, Rjb, vs30,U=U,SS=SS,RS=RS,NS=NS,Z1=None,intensity_measure='PGV',italy=True)
Sa03_bssa14,stdev_Sa03=bssa14(M, Rjb, vs30,U=U,SS=SS,RS=RS,NS=NS,Z1=None,intensity_measure='Sa0.3',italy=True)
Sa10_bssa14,stdev_Sa10=bssa14(M, Rjb, vs30,U=U,SS=SS,RS=RS,NS=NS,Z1=None,intensity_measure='Sa1.0',italy=True)
Sa30_bssa14,stdev_Sa30=bssa14(M, Rjb, vs30,U=U,SS=SS,RS=RS,NS=NS,Z1=None,intensity_measure='Sa3.0',italy=True)

#Stdevs
pga_plus=exp(log(pga_bssa14)+stdev_pga)
pga_minus=exp(log(pga_bssa14)-stdev_pga)
pgv_plus=exp(log(pgv_bssa14)+stdev_pgv)
pgv_minus=exp(log(pgv_bssa14)-stdev_pgv)
Sa03_plus=exp(log(Sa03_bssa14)+stdev_Sa03)
Sa03_minus=exp(log(Sa03_bssa14)-stdev_Sa03)
Sa10_plus=exp(log(Sa10_bssa14)+stdev_Sa10)
Sa10_minus=exp(log(Sa10_bssa14)-stdev_Sa10)
Sa30_plus=exp(log(Sa30_bssa14)+stdev_Sa30)
Sa30_minus=exp(log(Sa30_bssa14)-stdev_Sa30)


plt.figure(figsize=(8,5))

plt.subplot(231)
i=where(n_or_s>0)
plt.scatter(Rjb_obs[i],pga[i],c='#1E90FF',marker='+',s=50)
i=where(n_or_s<0)
plt.scatter(Rjb_obs[i],pga[i],c='#DC143C',marker='+',s=50)
plt.loglog(Rjb,pga_bssa14,lw=2,c='#606060')
plt.loglog(Ritaly,pga_italy,lw=2,c='#228B22')
plt.loglog(Rjb,pga_minus,lw=1,c='#606060')
plt.loglog(Rjb,pga_plus,lw=1,c='#606060')
plt.ylim([0.002,1.6])
plt.xlim([1,200])
plt.annotate('PGA (g)',xy=(2,0.005))
ax=plt.gca()
ax.xaxis.set_ticklabels([])



plt.subplot(232)
i=where(n_or_s>0)
plt.scatter(Rjb_obs[i],pgv[i],c='#1E90FF',marker='+',s=50)
i=where(n_or_s<0)
plt.scatter(Rjb_obs[i],pgv[i],c='#DC143C',marker='+',s=50)
plt.loglog(Ritaly,pgv_italy,lw=2,c='#228B22')
plt.loglog(Rjb,pgv_bssa14,lw=2,c='#606060')
plt.legend(['BSSA14','MA2011','North of hypo','South of hypo'],loc=3,frameon=False,bbox_to_anchor=[1.0, 0])
plt.loglog(Rjb,pgv_minus,lw=1,c='#606060')
plt.loglog(Rjb,pgv_plus,lw=1,c='#606060')
plt.ylim([0.3,90])
plt.xlim([1,200])
plt.annotate('PGV (cm/s)',xy=(2,0.6))
ax=plt.gca()
ax.xaxis.set_ticklabels([])



plt.subplot(234)
plt.loglog(Rjb,Sa03_bssa14,lw=2,c='#606060')
plt.loglog(Rjb,Sa03_minus,lw=1,c='#606060')
plt.loglog(Rjb,Sa03_plus,lw=1,c='#606060')
plt.loglog(Ritaly,Sa03_italy,lw=2,c='#228B22')
i=where(n_or_s>0)
plt.scatter(Rjb_obs[i],Sa03[i],c='#1E90FF',marker='+',s=50)
i=where(n_or_s<0)
plt.scatter(Rjb_obs[i],Sa03[i],c='#DC143C',marker='+',s=50)
plt.ylim([0.002,1.6])
plt.xlim([1,200])
plt.xlabel('Distance (km)')
plt.annotate('Sa 0.3s (g)',xy=(2,0.005))


plt.subplot(235)
plt.loglog(Rjb,Sa10_bssa14,lw=2,c='#606060')
plt.loglog(Rjb,Sa10_minus,lw=1,c='#606060')
plt.loglog(Rjb,Sa10_plus,lw=1,c='#606060')
plt.loglog(Ritaly,Sa10_italy,lw=2,c='#228B22')
i=where(n_or_s>0)
plt.scatter(Rjb_obs[i],Sa10[i],c='#1E90FF',marker='+',s=50)
i=where(n_or_s<0)
plt.scatter(Rjb_obs[i],Sa10[i],c='#DC143C',marker='+',s=50)
plt.ylim([0.002,1.6])
plt.xlim([1,200])
plt.xlabel('Distance (km)')
plt.annotate('Sa 1.0s (g)',xy=(2,0.005))


plt.subplot(236)
plt.loglog(Rjb,Sa30_bssa14,lw=2,c='#606060')
plt.loglog(Rjb,Sa30_minus,lw=1,c='#606060')
plt.loglog(Rjb,Sa30_plus,lw=1,c='#606060')
plt.loglog(Ritaly,Sa30_italy,lw=2,c='#228B22')
i=where(n_or_s>0)
plt.scatter(Rjb_obs[i],Sa30[i],c='#1E90FF',marker='+',s=50)
i=where(n_or_s<0)
plt.scatter(Rjb_obs[i],Sa30[i],c='#DC143C',marker='+',s=50)
plt.ylim([0.002,1.6])
plt.xlim([1,200])
plt.xlabel('Distance (km)')
plt.annotate('Sa 3.0s (g)',xy=(2,0.005))

plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.95,hspace=0.02)

plt.show()
