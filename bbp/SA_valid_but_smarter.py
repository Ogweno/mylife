from numpy import genfromtxt,zeros,r_,diff
from obspy import read
from glob import glob
from scipy.integrate import cumtrapz
from matplotlib import pyplot as plt
from numpy import reshape,array
from mudpy.ruptfunctions import rotatedResponseSpectrum

#Process fakequakes data for PGA and PGV

Nsources=100
path='/Users/dmelgar/FakeQuakes/M6_validation/output/waveforms/'
folders=glob(path+'M6*')
stations=['JPSB','LUTZ','MHCB','MILP','MONB']

bbfolders=glob('/Users/dmelgar/FakeQuakes/M6_validation/bbp_data/*')

SA_fq=zeros((Nsources,5))
SA_bbp=zeros((Nsources,5))

fig, axarr = plt.subplots(5, 5,figsize=(17,12))

f=array([0.1,0.3,1.0,3.0,10.0])

for ksta in range(len(stations)):
    
    print stations[ksta]

    for k in range(Nsources):
            
        an=read(folders[k]+'/'+stations[ksta]+'.bb.HNN.sac')
        ae=read(folders[k]+'/'+stations[ksta]+'.bb.HNE.sac')
        
        #GET SA
        SA=rotatedResponseSpectrum(ae[0].stats.delta,ae[0].data/9.81,an[0].data/9.81,f)
        SA_fq[k,:]=SA[0]['value']
        
            
        ##### now get BBP stuff
        bbfolders[k]
        root=bbfolders[k].split('/')[-1]
        sta=stations[ksta]
        a=genfromtxt(bbfolders[k]+'/'+root+'.'+sta+'.acc.bbp')

        #GET SA
        dt=diff(a[:,0])[0]
        SA=rotatedResponseSpectrum(0.025,a[:,2]/981,a[:,1]/981,f)
        SA_bbp[k,:]=SA[0]['value']
    
    
    #####    SA 0.1    ##### 
    ax=axarr[0,ksta] #current axis   
    
    x=range(Nsources)
    kf=0

    ax.plot([1,Nsources],[SA_fq[:,kf].mean(),SA_fq[:,kf].mean()],'r',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean(),SA_bbp[:,kf].mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('SA(f=0.1Hz)g')
        ax.legend(['FQ','BBP'],frameon=False)
    
    ax.plot(x,SA_fq[:,kf],'r',lw=0.5)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()+SA_fq[:,kf].std(),SA_fq[:,kf].mean()+SA_fq[:,kf].std()],'r',lw=2)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()-SA_fq[:,kf].std(),SA_fq[:,kf].mean()-SA_fq[:,kf].std()],'r',lw=2)
    
    ax.plot(x,SA_bbp[:,kf],'b',lw=0.5)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()+SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()+SA_bbp[:,kf].std()],'b',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()-SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()-SA_bbp[:,kf].std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,1.1*yl[1]])
    ax.set_title(stations[ksta])
    

    
        
    #####    SA 0.3     #####    
    
    kf=1
    ax=axarr[kf,ksta] #current axis 

    ax.plot([1,Nsources],[SA_fq[:,kf].mean(),SA_fq[:,kf].mean()],'r',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean(),SA_bbp[:,kf].mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('SA(f=0.3Hz)g')
    
    ax.plot(x,SA_fq[:,kf],'r',lw=0.5)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()+SA_fq[:,kf].std(),SA_fq[:,kf].mean()+SA_fq[:,kf].std()],'r',lw=2)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()-SA_fq[:,kf].std(),SA_fq[:,kf].mean()-SA_fq[:,kf].std()],'r',lw=2)
    
    ax.plot(x,SA_bbp[:,kf],'b',lw=0.5)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()+SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()+SA_bbp[:,kf].std()],'b',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()-SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()-SA_bbp[:,kf].std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,1.1*yl[1]])
    
    
    #####    SA 1.0     #####    
    
    kf=2
    ax=axarr[kf,ksta] #current axis 

    ax.plot([1,Nsources],[SA_fq[:,kf].mean(),SA_fq[:,kf].mean()],'r',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean(),SA_bbp[:,kf].mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('SA(f=1.0Hz)g')
    
    ax.plot(x,SA_fq[:,kf],'r',lw=0.5)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()+SA_fq[:,kf].std(),SA_fq[:,kf].mean()+SA_fq[:,kf].std()],'r',lw=2)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()-SA_fq[:,kf].std(),SA_fq[:,kf].mean()-SA_fq[:,kf].std()],'r',lw=2)
    
    ax.plot(x,SA_bbp[:,kf],'b',lw=0.5)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()+SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()+SA_bbp[:,kf].std()],'b',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()-SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()-SA_bbp[:,kf].std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,1.1*yl[1]])
    
    
    #####    SA 3.0    #####    
    
    kf=3
    ax=axarr[kf,ksta] #current axis 

    ax.plot([1,Nsources],[SA_fq[:,kf].mean(),SA_fq[:,kf].mean()],'r',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean(),SA_bbp[:,kf].mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('SA(f=3.0Hz)g')
    
    ax.plot(x,SA_fq[:,kf],'r',lw=0.5)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()+SA_fq[:,kf].std(),SA_fq[:,kf].mean()+SA_fq[:,kf].std()],'r',lw=2)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()-SA_fq[:,kf].std(),SA_fq[:,kf].mean()-SA_fq[:,kf].std()],'r',lw=2)
    
    ax.plot(x,SA_bbp[:,kf],'b',lw=0.5)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()+SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()+SA_bbp[:,kf].std()],'b',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()-SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()-SA_bbp[:,kf].std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,1.1*yl[1]])
 
        
    #####    SA 10.0    #####    
    
    kf=4
    ax=axarr[kf,ksta] #current axis 

    ax.plot([1,Nsources],[SA_fq[:,kf].mean(),SA_fq[:,kf].mean()],'r',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean(),SA_bbp[:,kf].mean()],'b',lw=2)
    if ksta==0:
        ax.set_ylabel('SA(f=10.0Hz)g')
    
    ax.plot(x,SA_fq[:,kf],'r',lw=0.5)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()+SA_fq[:,kf].std(),SA_fq[:,kf].mean()+SA_fq[:,kf].std()],'r',lw=2)
    ax.plot([1,Nsources],[SA_fq[:,kf].mean()-SA_fq[:,kf].std(),SA_fq[:,kf].mean()-SA_fq[:,kf].std()],'r',lw=2)
    
    ax.plot(x,SA_bbp[:,kf],'b',lw=0.5)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()+SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()+SA_bbp[:,kf].std()],'b',lw=2)
    ax.plot([1,Nsources],[SA_bbp[:,kf].mean()-SA_bbp[:,kf].std(),SA_bbp[:,kf].mean()-SA_bbp[:,kf].std()],'b',lw=2)
    
    yl=ax.get_ylim()
    ax.set_ylim([0,1.1*yl[1]])  
    ax.set_xlabel('Realization')     

plt.subplots_adjust(left=0.06,right=0.98,top=0.95,bottom=0.1)


plt.show()