from numpy import genfromtxt,zeros,r_,diff
from obspy import read
from glob import glob
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt
from numpy import reshape

#Process fakequakes data for PGA and PGV

Nsources=100
path='/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/'
folders=glob(path+'M6*')
stations=['JPSB','LUTZ','MHCB','MILP','MONB']

bbfolders=glob('/Users/dmelgar/FakeQuakes/M6_validation/bbp_data/*')

pga_fq=zeros(Nsources)
pgv_fq=zeros(Nsources)
pga_bbp=zeros(Nsources)
pgv_bbp=zeros(Nsources)


fig, axarr = plt.subplots(2, 5,figsize=(16,9))

for ksta in range(len(stations)):
    
    print stations[ksta]

    for k in range(Nsources):
            
        an=read(folders[k]+'/'+stations[ksta]+'.bb.HNN.sac')
        ae=read(folders[k]+'/'+stations[ksta]+'.bb.HNE.sac')
        
        vn=an.copy()
        ve=ae.copy()
        
        vn[0].data=cumtrapz(an[0].data,an[0].times(),initial=0)
        ve[0].data=cumtrapz(ae[0].data,ae[0].times(),initial=0)
        
        #Get peak values
        pga_fq[k]=max(r_[abs(an[0].data),abs(ae[0].data)])
        pgv_fq[k]=max(r_[abs(vn[0].data),abs(ve[0].data)])
        
            
        ##### now get BBP stuff
        bbfolders[k]
        root=bbfolders[k].split('/')[-1]
        sta=stations[ksta]
        a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')
        v=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.vel.bbp')
        
        pga_bbp[k]=abs(a[:,1:3]).max()/100.
        pgv_bbp[k]=abs(v[:,1:3]).max()/100.
    
    #Make plot
    ax=axarr[0,ksta] #current axis   
    
    x=range(Nsources)

    ax.plot([1,Nsources],[pga_fq.mean(),pga_fq.mean()],'r',lw=2)
    ax.plot([1,Nsources],[pga_bbp.mean(),pga_bbp.mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('PGA (m/s/s)')
        ax.legend(['FQ','BBP'],frameon=False)
    
    ax.plot(x,pga_fq,'r',lw=0.5)
    ax.plot([1,Nsources],[pga_fq.mean()+pga_fq.std(),pga_fq.mean()+pga_fq.std()],'r',lw=2)
    ax.plot([1,Nsources],[pga_fq.mean()-pga_fq.std(),pga_fq.mean()-pga_fq.std()],'r',lw=2)
    
    ax.plot(x,pga_bbp,'b',lw=0.5)
    ax.plot([1,Nsources],[pga_bbp.mean()+pga_bbp.std(),pga_bbp.mean()+pga_bbp.std()],'b',lw=2)
    ax.plot([1,Nsources],[pga_bbp.mean()-pga_bbp.std(),pga_bbp.mean()-pga_bbp.std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,yl[1]])
    ax.set_title(stations[ksta])
    

    
        
    #####    PGV      #####    
    
    ax=axarr[1,ksta] #current axis   
    
    ax.plot(x,pgv_fq,'r',lw=0.5)
    ax.plot([1,Nsources],[pgv_fq.mean(),pgv_fq.mean()],'r',lw=2)
    ax.plot([1,Nsources],[pgv_fq.mean()+pgv_fq.std(),pgv_fq.mean()+pgv_fq.std()],'r',lw=2)
    ax.plot([1,Nsources],[pgv_fq.mean()-pgv_fq.std(),pgv_fq.mean()-pgv_fq.std()],'r',lw=2)
    
    ax.plot(x,pgv_bbp,'b',lw=0.5)
    ax.plot([1,Nsources],[pgv_bbp.mean(),pgv_bbp.mean()],'b',lw=2)
    ax.plot([1,Nsources],[pgv_bbp.mean()+pgv_bbp.std(),pgv_bbp.mean()+pgv_bbp.std()],'b',lw=2)
    ax.plot([1,Nsources],[pgv_bbp.mean()-pgv_bbp.std(),pgv_bbp.mean()-pgv_bbp.std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,yl[1]])
    ax.set_xlabel('Realization')
    
    if ksta==0:
        ax.set_ylabel('PGV (m/s)')
        

plt.subplots_adjust(left=0.06,right=0.98,top=0.95,bottom=0.1)


plt.show()