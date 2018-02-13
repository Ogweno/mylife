from numpy import genfromtxt,r_,zeros,unique,where,diff
from matplotlib import pyplot as plt

variable='Mxx'
f0=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process0')
f1=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process1')
f2=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process2')
f3=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process3')
f4=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process4')
f5=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process5')
f6=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process6')
f7=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/tmp_'+variable+'_process7')

f=r_[f0,f1,f2,f3,f4,f5,f6,f7]
s=unique(f[:,1])
count=zeros(len(s))
for k in range(len(s)):
    print k
    i=where(f[:,1]==s[k])[0]
    count[k]=len(i)
    
    
plt.figure()
plt.scatter(s,count)
plt.plot(diff(s))