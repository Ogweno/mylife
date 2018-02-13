from glob import glob
from numpy import genfromtxt,where,log10,r_,ones,array,logspace,argmin,arange
import bbptools
from scipy.integrate import cumtrapz
from pyproj import Geod
from matplotlib import pyplot as plt
from mudpy.forward import highpass

magnitude='m6.5'
M=6.5
yl=[1e-3,5e1]
path='/Users/dmelgar/code/BBP/bbp/bbp_data/finished/large_eew_test/'
folders=glob(path+'*'+magnitude+'*')
filter=False
corner=1./13

#Initalize
observed_Pd=[]
predicted_Pd=[]
M=array([])
R=array([])

#Loop over sim_list
for ksim in folders:
    print ksim
    stations=genfromtxt(glob(ksim+'/param_files/*.stl')[0],usecols=2,dtype='S')
    lonlat=genfromtxt(glob(ksim+'/param_files/*.stl')[0],usecols=[0,1])
    
    #Get ID#
    files=glob(ksim+'/*.bbp')
    ID=files[0].split('/')[-1].split('.')[0]
    
    for kstation in range(len(stations)):
        
        current_station=stations[kstation]
        
        v=genfromtxt(ksim+'/'+ID+'.'+current_station+'.vel.bbp')
        t=v[:,0]
        vz=v[:,3]
        #Get displacement
        dz=cumtrapz(vz,t,initial=0)
        
        dt=t[1]-t[0]
        if filter==True:
            dfil=highpass(dz,corner,1/dt,4)
            dz=dfil.copy()
        
        #Read meta data get magnitude and hypocenter
        if kstation==0:
            #Get hypocenter
            srf_file=glob(ksim+'/*.srf')[0]
            xyz,slip,tinit,stf_all=bbptools.read_srf(srf_file)
            i=argmin(tinit)
            hypocenter=xyz[i,:]
            print hypocenter

            #get magnitude
            metadata=glob(ksim+'/param_files/*.src')[0]          
            f=open(metadata,'r')
            line=f.readline()
            mag=float(line.split()[-1])
            f.close()
            M=r_[M,mag*ones(len(stations))]
        
        #Get P,S arrivals
        ptime,stime=bbptools.arrivals(hypocenter,lonlat[kstation,0],lonlat[kstation,1])
        
        #Get  observed Pd
        ptime=ptime
        i=where((t>=ptime) & (t<=ptime+3))[0]
        Pd=max(abs(dz[i]))
        observed_Pd.append(Pd)
        
        #Get predicted Pd
        
        #Station to hypo distance
        g=Geod(ellps='WGS84')
        az,baz,dist=g.inv(hypocenter[0],hypocenter[1],lonlat[kstation,0],lonlat[kstation,1])
        #Convert distance from m to degrees and km
        dist_km=dist/1000
        R=r_[R,dist_km]
            
        
        #Finally, actually get Pd
        Pd=10**((mag-5.39-1.38*log10(dist_km))/1.23)      
        predicted_Pd.append(Pd)  
        
observed_Pd=array(observed_Pd)
predicted_Pd=array(predicted_Pd)
        
#Make plot
plt.figure(figsize=(5,5))    
k=0
i=arange(46*k,46*(k+1))
plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DC143C')
k=1
i=arange(46*k,46*(k+1))
plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#32CD32')
k=2
i=arange(46*k,46*(k+1))
plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#0000CD')
k=3
i=arange(46*k,46*(k+1))
plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#DAA520')
k=4
i=arange(46*k,46*(k+1))
plt.scatter(R[i],observed_Pd[i],lw=0.5,s=40,c='#9932CC')


ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')



ax.set_ylim(yl)
ax.set_xlim([1e0,1e2])
plt.xlabel('Distance (km)')
plt.ylabel('Pd (cm)')

plt.legend(['5%','25%','50%','75%','95%'],loc=3)
plt.title(magnitude.upper())

#Reference lines
Rref=logspace(0,2)
Pd=10**((M[0]-5.39-1.38*log10(Rref))/1.23)
Pd_plus=10**((M[0]-0.3-5.39-1.38*log10(Rref))/1.23)
Pd_minus=10**((M[0]+0.3-5.39-1.38*log10(Rref))/1.23)

plt.plot(Rref,Pd,'--',lw=2,c='#505050')
plt.plot(Rref,Pd_plus,'--',lw=2,c='#505050')
plt.plot(Rref,Pd_minus,'--',lw=2,c='#505050')

plt.subplots_adjust(left=0.15,right=0.97,top=0.94,bottom=0.11)

plt.show()