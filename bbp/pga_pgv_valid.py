from numpy import genfromtxt,zeros,r_,diff
from obspy import read
from glob import glob
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt

#Process fakequakes data for PGA and PGV

path='/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/'
folders=glob(path+'M6*')
stations=['JPSB','LUTZ','MHCB','MILP','MONB']

bbfolders=glob('/Users/dmelgar/FakeQuakes/M6_validation/bbp_results/outdata/*')

jpsb_pga=zeros((20,3))
jpsb_pgv=zeros((20,3))
lutz_pga=zeros((20,3))
lutz_pgv=zeros((20,3))
mhcb_pga=zeros((20,3))
mhcb_pgv=zeros((20,3))
milp_pga=zeros((20,3))
milp_pgv=zeros((20,3))
monp_pga=zeros((20,3))
monp_pgv=zeros((20,3))

jpsb_pga_bbp=zeros((20,3))
jpsb_pgv_bbp=zeros((20,3))
lutz_pga_bbp=zeros((20,3))
lutz_pgv_bbp=zeros((20,3))
mhcb_pga_bbp=zeros((20,3))
mhcb_pgv_bbp=zeros((20,3))
milp_pga_bbp=zeros((20,3))
milp_pgv_bbp=zeros((20,3))
monp_pga_bbp=zeros((20,3))
monp_pgv_bbp=zeros((20,3))

for k in range(20):
    
    an=read(folders[k]+'/'+stations[0]+'.bb.HNN.sac')
    ae=read(folders[k]+'/'+stations[0]+'.bb.HNE.sac')
    az=read(folders[k]+'/'+stations[0]+'.bb.HNZ.sac')
    
    vn=an.copy()
    ve=ae.copy()
    vz=az.copy()
    
    vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    #Get peak values
    jpsb_pga[k,0]=max(abs(an[0].data))
    jpsb_pga[k,1]=max(abs(ae[0].data))
    jpsb_pga[k,2]=max(abs(az[0].data))
    jpsb_pgv[k,0]=max(abs(vn[0].data))
    jpsb_pgv[k,1]=max(abs(ve[0].data))
    jpsb_pgv[k,2]=max(abs(vz[0].data))
    
    ######    
    an=read(folders[k]+'/'+stations[1]+'.bb.HNN.sac')
    ae=read(folders[k]+'/'+stations[1]+'.bb.HNE.sac')
    az=read(folders[k]+'/'+stations[1]+'.bb.HNZ.sac')
    
    vn=an.copy()
    ve=ae.copy()
    vz=az.copy()
    
    vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    #Get peak values
    lutz_pga[k,0]=max(abs(an[0].data))
    lutz_pga[k,1]=max(abs(ae[0].data))
    lutz_pga[k,2]=max(abs(az[0].data))
    lutz_pgv[k,0]=max(abs(vn[0].data))
    lutz_pgv[k,1]=max(abs(ve[0].data))
    lutz_pgv[k,2]=max(abs(vz[0].data))
    
    ######    
    an=read(folders[k]+'/'+stations[2]+'.bb.HNN.sac')
    ae=read(folders[k]+'/'+stations[2]+'.bb.HNE.sac')
    az=read(folders[k]+'/'+stations[2]+'.bb.HNZ.sac')
    
    vn=an.copy()
    ve=ae.copy()
    vz=az.copy()
    
    vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    #Get peak values
    mhcb_pga[k,0]=max(abs(an[0].data))
    mhcb_pga[k,1]=max(abs(ae[0].data))
    mhcb_pga[k,2]=max(abs(az[0].data))
    mhcb_pgv[k,0]=max(abs(vn[0].data))
    mhcb_pgv[k,1]=max(abs(ve[0].data))
    mhcb_pgv[k,2]=max(abs(vz[0].data))
    
    ######    
    an=read(folders[k]+'/'+stations[3]+'.bb.HNN.sac')
    ae=read(folders[k]+'/'+stations[3]+'.bb.HNE.sac')
    az=read(folders[k]+'/'+stations[3]+'.bb.HNZ.sac')
    
    vn=an.copy()
    ve=ae.copy()
    vz=az.copy()
    
    vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    #Get peak values
    milp_pga[k,0]=max(abs(an[0].data))
    milp_pga[k,1]=max(abs(ae[0].data))
    milp_pga[k,2]=max(abs(az[0].data))
    milp_pgv[k,0]=max(abs(vn[0].data))
    milp_pgv[k,1]=max(abs(ve[0].data))
    milp_pgv[k,2]=max(abs(vz[0].data))
    
    ######    
    an=read(folders[k]+'/'+stations[4]+'.bb.HNN.sac')
    ae=read(folders[k]+'/'+stations[4]+'.bb.HNE.sac')
    az=read(folders[k]+'/'+stations[4]+'.bb.HNZ.sac')
    
    vn=an.copy()
    ve=ae.copy()
    vz=az.copy()
    
    vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
    ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
    vz[0].data=cumtrapz(az[0].data,az[0].times(),initial=0)
    
    #Get peak values
    monp_pga[k,0]=max(abs(an[0].data))
    monp_pga[k,1]=max(abs(ae[0].data))
    monp_pga[k,2]=max(abs(az[0].data))
    monp_pgv[k,0]=max(abs(vn[0].data))
    monp_pgv[k,1]=max(abs(ve[0].data))
    monp_pgv[k,2]=max(abs(vz[0].data))
    
    
    
    ##### now get BBP stuff
    bbfolders[k]
    root=bbfolders[k].split('/')[-1]
    sta=stations[0]
    a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
    v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
    
    jpsb_pga_bbp[k,0]=max(abs(a[:,1]))/100.
    jpsb_pga_bbp[k,1]=max(abs(a[:,2]))/100.
    jpsb_pga_bbp[k,2]=max(abs(a[:,3]))/100.
    jpsb_pgv_bbp[k,0]=max(abs(v[:,1]))/100.
    jpsb_pgv_bbp[k,1]=max(abs(v[:,2]))/100.
    jpsb_pgv_bbp[k,2]=max(abs(v[:,3]))/100.
    
    ####
    root=bbfolders[k].split('/')[-1]
    sta=stations[1]
    a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
    v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
    
    lutz_pga_bbp[k,0]=max(abs(a[:,1]))/100.
    lutz_pga_bbp[k,1]=max(abs(a[:,2]))/100.
    lutz_pga_bbp[k,2]=max(abs(a[:,3]))/100.
    lutz_pgv_bbp[k,0]=max(abs(v[:,1]))/100.
    lutz_pgv_bbp[k,1]=max(abs(v[:,2]))/100.
    lutz_pgv_bbp[k,2]=max(abs(v[:,3]))/100.
    
    
    ####
    root=bbfolders[k].split('/')[-1]
    sta=stations[2]
    a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
    v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
    
    mhcb_pga_bbp[k,0]=max(abs(a[:,1]))/100.
    mhcb_pga_bbp[k,1]=max(abs(a[:,2]))/100.
    mhcb_pga_bbp[k,2]=max(abs(a[:,3]))/100.
    mhcb_pgv_bbp[k,0]=max(abs(v[:,1]))/100.
    mhcb_pgv_bbp[k,1]=max(abs(v[:,2]))/100.
    mhcb_pgv_bbp[k,2]=max(abs(v[:,3]))/100.
    
    
    ####
    root=bbfolders[k].split('/')[-1]
    sta=stations[3]
    a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
    v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
    
    milp_pga_bbp[k,0]=max(abs(a[:,1]))/100.
    milp_pga_bbp[k,1]=max(abs(a[:,2]))/100.
    milp_pga_bbp[k,2]=max(abs(a[:,3]))/100.
    milp_pgv_bbp[k,0]=max(abs(v[:,1]))/100.
    milp_pgv_bbp[k,1]=max(abs(v[:,2]))/100.
    milp_pgv_bbp[k,2]=max(abs(v[:,3]))/100.
    
    
    ####
    root=bbfolders[k].split('/')[-1]
    sta=stations[4]
    a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
    v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
    
    monp_pga_bbp[k,0]=max(abs(a[:,1]))/100.
    monp_pga_bbp[k,1]=max(abs(a[:,2]))/100.
    monp_pga_bbp[k,2]=max(abs(a[:,3]))/100.
    monp_pgv_bbp[k,0]=max(abs(v[:,1]))/100.
    monp_pgv_bbp[k,1]=max(abs(v[:,2]))/100.
    monp_pgv_bbp[k,2]=max(abs(v[:,3]))/100.
    
    
#Get maximum
    
x=range(20)
plt.figure()
plt.subplot(251)
plt.plot(x,jpsb_pga[:,0],'r')
plt.plot(x,jpsb_pga_bbp[:,0],'--r')
plt.plot(x,jpsb_pga[:,1],'g')
plt.plot(x,jpsb_pga_bbp[:,1],'--g')
plt.plot(x,jpsb_pga[:,2],'b')
plt.plot(x,jpsb_pga_bbp[:,2],'--b')
plt.ylabel('PGA m/s/s')
plt.title('JPSB')

plt.subplot(252)
plt.plot(x,lutz_pga[:,0],'r')
plt.plot(x,lutz_pga_bbp[:,0],'--r')
plt.plot(x,lutz_pga[:,1],'g')
plt.plot(x,lutz_pga_bbp[:,1],'--g')
plt.plot(x,lutz_pga[:,2],'b')
plt.plot(x,lutz_pga_bbp[:,2],'--b')
plt.title('LUTZ')

plt.subplot(253)
plt.plot(x,mhcb_pga[:,0],'r')
plt.plot(x,mhcb_pga_bbp[:,0],'--r')
plt.plot(x,mhcb_pga[:,1],'g')
plt.plot(x,mhcb_pga_bbp[:,1],'--g')
plt.plot(x,mhcb_pga[:,2],'b')
plt.plot(x,mhcb_pga_bbp[:,2],'--b')
plt.title('MHCB')

plt.subplot(254)
plt.plot(x,milp_pga[:,0],'r')
plt.plot(x,milp_pga_bbp[:,0],'--r')
plt.plot(x,milp_pga[:,1],'g')
plt.plot(x,milp_pga_bbp[:,1],'--g')
plt.plot(x,milp_pga[:,2],'b')
plt.plot(x,milp_pga_bbp[:,2],'--b')
plt.title('MILP')

plt.subplot(255)
plt.plot(x,monp_pga[:,0],'r')
plt.plot(x,monp_pga_bbp[:,0],'--r')
plt.plot(x,monp_pga[:,1],'g')
plt.plot(x,monp_pga_bbp[:,1],'--g')
plt.plot(x,monp_pga[:,2],'b')
plt.plot(x,monp_pga_bbp[:,2],'--b')
plt.title('MONP')



plt.subplot(256)
plt.plot(x,jpsb_pgv[:,0],'r')
plt.plot(x,jpsb_pgv_bbp[:,0],'--r')
plt.plot(x,jpsb_pgv[:,1],'g')
plt.plot(x,jpsb_pgv_bbp[:,1],'--g')
plt.plot(x,jpsb_pgv[:,2],'b')
plt.plot(x,jpsb_pgv_bbp[:,2],'--b')
plt.ylabel('PGV m/s')
plt.title('JPSB')

plt.subplot(257)
plt.plot(x,lutz_pgv[:,0],'r')
plt.plot(x,lutz_pgv_bbp[:,0],'--r')
plt.plot(x,lutz_pgv[:,1],'g')
plt.plot(x,lutz_pgv_bbp[:,1],'--g')
plt.plot(x,lutz_pgv[:,2],'b')
plt.plot(x,lutz_pgv_bbp[:,2],'--b')
plt.title('LUTZ')

plt.subplot(258)
plt.plot(x,mhcb_pgv[:,0],'r')
plt.plot(x,mhcb_pgv_bbp[:,0],'--r')
plt.plot(x,mhcb_pgv[:,1],'g')
plt.plot(x,mhcb_pgv_bbp[:,1],'--g')
plt.plot(x,mhcb_pgv[:,2],'b')
plt.plot(x,mhcb_pgv_bbp[:,2],'--b')
plt.title('MHCB')

plt.subplot(259)
plt.plot(x,milp_pgv[:,0],'r')
plt.plot(x,milp_pgv_bbp[:,0],'--r')
plt.plot(x,milp_pgv[:,1],'g')
plt.plot(x,milp_pgv_bbp[:,1],'--g')
plt.plot(x,milp_pgv[:,2],'b')
plt.plot(x,milp_pgv_bbp[:,2],'--b')
plt.title('MILP')

plt.subplot(2,5,10)
plt.plot(x,monp_pgv[:,0],'r')
plt.plot(x,monp_pgv_bbp[:,0],'--r')
plt.plot(x,monp_pgv[:,1],'g')
plt.plot(x,monp_pgv_bbp[:,1],'--g')
plt.plot(x,monp_pgv[:,2],'b')
plt.plot(x,monp_pgv_bbp[:,2],'--b')
plt.title('MONP')



x=range(20)
plt.figure()
plt.subplot(251)
plt.plot(x,jpsb_pga[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(jpsb_pga[:,0:2].max(axis=1)),mean(jpsb_pga[:,0:2].max(axis=1)),],'r')
plt.plot(x,jpsb_pga_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(jpsb_pga_bbp[:,0:2].max(axis=1)),mean(jpsb_pga_bbp[:,0:2].max(axis=1)),],'b')
plt.ylabel('PGA m/s/s')
plt.title('JPSB')

plt.subplot(252)
plt.plot(x,lutz_pga[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(lutz_pga[:,0:2].max(axis=1)),mean(lutz_pga[:,0:2].max(axis=1)),],'r')
plt.plot(x,lutz_pga_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(lutz_pga_bbp[:,0:2].max(axis=1)),mean(lutz_pga_bbp[:,0:2].max(axis=1)),],'b')
plt.title('LUTZ')

plt.subplot(253)
plt.plot(x,mhcb_pga[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(mhcb_pga[:,0:2].max(axis=1)),mean(mhcb_pga[:,0:2].max(axis=1)),],'r')
plt.plot(x,mhcb_pga_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(mhcb_pga_bbp[:,0:2].max(axis=1)),mean(mhcb_pga_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MHCB')

plt.subplot(254)
plt.plot(x,milp_pga[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(milp_pga[:,0:2].max(axis=1)),mean(milp_pga[:,0:2].max(axis=1)),],'r')
plt.plot(x,milp_pga_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(milp_pga_bbp[:,0:2].max(axis=1)),mean(milp_pga_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MILP')

plt.subplot(255)
plt.plot(x,monp_pga[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(monp_pga[:,0:2].max(axis=1)),mean(monp_pga[:,0:2].max(axis=1)),],'r')
plt.plot(x,monp_pga_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(monp_pga_bbp[:,0:2].max(axis=1)),mean(monp_pga_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MONP')


plt.subplot(256)
plt.plot(x,jpsb_pgv[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(jpsb_pgv[:,0:2].max(axis=1)),mean(jpsb_pgv[:,0:2].max(axis=1)),],'r')
plt.plot(x,jpsb_pgv_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(jpsb_pgv_bbp[:,0:2].max(axis=1)),mean(jpsb_pgv_bbp[:,0:2].max(axis=1)),],'b')
plt.ylabel('PGV m/s/s')
plt.title('JPSB')

plt.subplot(257)
plt.plot(x,lutz_pgv[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(lutz_pgv[:,0:2].max(axis=1)),mean(lutz_pgv[:,0:2].max(axis=1)),],'r')
plt.plot(x,lutz_pgv_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(lutz_pgv_bbp[:,0:2].max(axis=1)),mean(lutz_pgv_bbp[:,0:2].max(axis=1)),],'b')
plt.title('LUTZ')

plt.subplot(258)
plt.plot(x,mhcb_pgv[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(mhcb_pgv[:,0:2].max(axis=1)),mean(mhcb_pgv[:,0:2].max(axis=1)),],'r')
plt.plot(x,mhcb_pgv_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(mhcb_pgv_bbp[:,0:2].max(axis=1)),mean(mhcb_pgv_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MHCB')

plt.subplot(259)
plt.plot(x,milp_pgv[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(milp_pgv[:,0:2].max(axis=1)),mean(milp_pgv[:,0:2].max(axis=1)),],'r')
plt.plot(x,milp_pgv_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(milp_pgv_bbp[:,0:2].max(axis=1)),mean(milp_pgv_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MILP')

plt.subplot(2,5,10)
plt.plot(x,monp_pgv[:,0:2].max(axis=1),'r')
plt.plot([1,20],[mean(monp_pgv[:,0:2].max(axis=1)),mean(monp_pgv[:,0:2].max(axis=1)),],'r')
plt.plot(x,monp_pgv_bbp[:,0:2].max(axis=1),'b')
plt.plot([1,20],[mean(monp_pgv_bbp[:,0:2].max(axis=1)),mean(monp_pgv_bbp[:,0:2].max(axis=1)),],'b')
plt.title('MONP')



plt.show()