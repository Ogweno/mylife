from matplotlib import pyplot as plt
from numpy import genfromtxt,where,argmin,ones,array,r_,arange,zeros,c_
from pyproj import Geod
from glob import glob
from scipy.interpolate import interp1d
from string import rjust

path_out='/Users/dmelgar/Chiapas2017/plots/tsunam_profile/'
eta_max=2.2
bathy=genfromtxt(u'/Users/dmelgar/Chiapas2017/misc/tsun_bathy_profile.txt')

# gauge files
path='/Users/dmelgar/Tsunamis/tehuantepec_profile/_output/'
files=glob(path+'gauge01*')


#Find trench
i=argmin(bathy[:,2])
trench=bathy[i,0:2]

#Get distances from bathy  to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(bathy[:,0],bathy[:,1],trench[0]*ones(len(bathy)),trench[1]*ones(len(bathy)))
dist=dist/1000
i=where(bathy[:,1]<trench[1])
dist[i]=-dist[i]
dist_bathy=dist
bathy=bathy[:,2]/1000

#time parameters
tout=arange(0,3600*24,60)
Ntimes=len(tout)
Ngauges=len(files)



#Read gauges interpolate an make output matrix
lon=zeros(len(files))
lat=zeros(len(files))
count=0
for k in range(len(files)):
    
    if k%20==0:
        print k
    
    f=open(files[k])
    line=f.readline()
    f.close()
    
    #read full gauge
    g=genfromtxt(files[k])
    g0=g[0,5]
    t=g[:,1]
    eta=g[:,5]
    
    #resample to correct times
    f=interp1d(t, eta, kind='linear')
    eta_interp=f(tout)
    
    #place in matrix
    if eta_interp.max()<eta_max:
        if count==0:
            eta_out=eta_interp
            lon=array(float(line.split()[4]))
            lat=array(float(line.split()[5]))
            
            count+=1
        else:
            eta_out=c_[eta_out,eta_interp]
            
            lon=c_[lon,float(line.split()[4])]
            lat=c_[lat,float(line.split()[5])]
    

#Get distances from tsunami_profile to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(lon[0,:],lat[0,:],trench[0]*ones(len(lon[0,:])),trench[1]*ones(len(lat[0,:])))
dist=dist/1000
i=where(lat[0,:]<trench[1])[0]
dist[i]=-dist[i]
dist_tsun=dist


#now for every time step plot profile
multiplier=4.5
offset=2
for kt in range(len(tout)):
    profile=eta_out[kt,:]
    print kt
    
    plt.figure(figsize=(16,6))
    plt.plot(dist_tsun,profile*multiplier+offset,lw=2,c='b')
    plt.plot(dist_bathy,bathy,lw=2,c='k')
    plt.xlabel('Distance from trench')
    plt.ylabel('Seafloor depth (m)')
    plt.plot([62.2,62.2],[-6.5,11],'--',lw=0.8)
    plt.plot([172,172],[-6.5,11],'--',lw=0.8)
    plt.xlim([-100,180])
    plt.ylim([-6.5,11])
    tmins=tout[kt]/60.
    plt.annotate(xy=(-90,8.5),s='t = %dmins' % (tmins),fontsize=16)
    
    plt.subplots_adjust(bottom=0.2)
    num=rjust(str(kt),4,'0')
    plt.savefig(path_out+num+'.png')
    plt.close()