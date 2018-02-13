from numpy import genfromtxt,where,arange,zeros,c_,r_,diff,nan,ones,savetxt,isnan,median
from glob import glob
from pyproj import Geod
from matplotlib import pyplot as plt

profiles=glob('/Users/dmelgar/DART_analysis/bathy_profiles/*dart.txt')
fout='/Users/dmelgar/DART_analysis/first_wave_to_dart.txt'

p=Geod(ellps='WGS84')

out=zeros((len(profiles),3))
for k in range(len(profiles)):
    print profiles[k]
    prof=genfromtxt(profiles[k])
    i=where(prof[:,2]<0)[0]
    
    if len(i)>2:
        prof=prof[i,:]
        az,baz,dx=p.inv(prof[0,0],prof[0,1],prof[1,0],prof[1,1])
        az,baz,dmax=p.inv(prof[0,0],prof[0,1],prof[-1,0],prof[-1,1])
        
        vel=(9.81*abs(prof[:,2]))**0.5 #in /m/s
        t=dx*ones(len(vel))/vel
        ttotal=t.sum()
        out[k,0:2]=prof[-1,0:2]
        out[k,2]=ttotal/60
    else:
        out[k,0:2]=prof[-1,1:3]
        out[k,2]=nan
        
savetxt(fout,out,fmt='%.4f')



bins=arange(0,240,10)
i=where(isnan(out[:,2])==False)[0]

plt.figure(figsize=(4,4))
plt.hist(out[i,2],bins=bins)
plt.plot([30,30],[0,160],'--')
#plt.plot([median(out[i,2]),median(out[i,2])],[0,160],'--')
plt.plot([median(out[i,2])/2,median(out[i,2])/2],[0,160],'--')
plt.subplots_adjust(bottom=0.13,left=0.17)
plt.xlabel('First arrival (mins)')
plt.ylabel('Count')
plt.ylim([0,160])

plt.show()