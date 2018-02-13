from matplotlib import pyplot as plt
import gmpe_tools
import bbptools
from numpy import logspace,zeros,ones,exp,log,log10,array,genfromtxt
from pyproj import Geod

path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/eew_test/'
sims=['45','50','55','60','65','70']
mags=['m4.5','m5.0','m5.5','m6.0','m6.5','m7.0']
xl=[1,50]

### PGA
plt.figure(figsize=(13,8))


#######      M4.5     ##########
k=0
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=4.5
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[0.01,10.0]



plt.subplot(231)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#DC143C')
plt.plot(R_pred,pred_pgd,lw=2,c='#DC143C')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#DC143C')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#DC143C')
plt.xlim(xl)
plt.ylim(yl)
plt.ylabel('PGD cm')
plt.legend(['M4.5'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M5.0     ##########
k=1
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=5.0
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[0.01,30.0]



plt.subplot(232)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#32CD32')
plt.plot(R_pred,pred_pgd,lw=2,c='#32CD32')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#32CD32')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#32CD32')
plt.xlim(xl)
plt.ylim(yl)
plt.legend(['M5.0'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M5.5     ##########
k=2
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=5.5
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[0.1,50.0]



plt.subplot(233)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#0000CD')
plt.plot(R_pred,pred_pgd,lw=2,c='#0000CD')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#0000CD')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#0000CD')
plt.xlim(xl)
plt.ylim(yl)
plt.legend(['M5.5'],loc=3,frameon=False)

########       END SUBPLOT     ########


#######      M6.0     ##########
k=3
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=6.0
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[0.3,120]



plt.subplot(234)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#DAA520')
plt.plot(R_pred,pred_pgd,lw=2,c='#DAA520')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#DAA520')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#DAA520')
plt.xlim(xl)
plt.ylim(yl)
plt.ylabel('PGD cm')
plt.xlabel('Hypo dist. (km)')
plt.legend(['M6.0'],loc=3,frameon=False)

########       END SUBPLOT     ########



#######      M6.5     ##########
k=4
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=6.5
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[2,1000.0]



plt.subplot(235)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#9932CC')
plt.plot(R_pred,pred_pgd,lw=2,c='#9932CC')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#9932CC')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#9932CC')
plt.xlim(xl)
plt.ylim(yl)
plt.xlabel('Hypo dist. (km)')
plt.legend(['M6.5'],loc=3,frameon=False)

########       END SUBPLOT     ########



#######      M7.0    ##########
k=5
sim_path=path+sims[k]
lonlat=genfromtxt(sim_path+'/param_files/eew_test.stl',usecols=[0,1])
stlon=lonlat[:,0]
stlat=lonlat[:,1]

#Get event hypocenter
srf_file=sim_path+'/template_'+mags[k]+'.srf'
xyz,slip,tinit,stf_all,rise_time,hypocenter=bbptools.read_srf(srf_file)

#Get station to event distances
evlon=ones(len(stlon))*hypocenter[0]
evlat=ones(len(stlat))*hypocenter[1]
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(stlon,stlat,evlon,evlat)
dist=dist/1000

#get predicted data
mag=7.0
delta=0.5
R_pred=logspace(-1,2)
pred_pgd=10**(-4.434+1.047*mag-0.138*mag*log10(R_pred))
pred_pgd_minus=10**(-4.434+1.047*mag-delta-0.138*(mag-delta)*log10(R_pred))
pred_pgd_plus=10**(-4.434+1.047*mag+delta-0.138*(mag+delta)*log10(R_pred))

#get simulation data
pga,pgv,pgd,Rjb,M=bbptools.get_simulation_intensity(sim_path)
pgdlims=log10(pgd)
pgd_range=pgdlims.max()-pgdlims.min()
yl=[6.0,4500.0]



plt.subplot(236)
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.scatter(dist,pgd,lw=0.5,s=40,c='#202020')
plt.plot(R_pred,pred_pgd,lw=2,c='#202020')
plt.plot(R_pred,pred_pgd_plus,'--',lw=2,c='#202020')
plt.plot(R_pred,pred_pgd_minus,'--',lw=2,c='#202020')
plt.xlim(xl)
plt.ylim(yl)
plt.xlabel('Hypo dist. (km)')
plt.legend(['M7.0'],loc=3,frameon=False)

########       END SUBPLOT     ########

plt.subplots_adjust(left=0.09,right=0.96,bottom=0.09,top=0.96,wspace=0.12,hspace=0.09)

plt.show()