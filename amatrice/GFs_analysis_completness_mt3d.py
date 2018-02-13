from numpy import genfromtxt,r_,zeros,unique,where,diff
from matplotlib import pyplot as plt

variable='Mxz'
f=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_'+variable+'.los')
#f=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_unfinished/_'+variable+'.unfinished3.los')
#f=genfromtxt(u'/Users/dmelgar/mt3d/Italy_modest_volume/GFs/static/_'+variable+'.los')

s=unique(f[:,1])
count=zeros(len(s))
for k in range(11600,len(s)):
    print k
    i=where(f[:,1]==s[k])[0]
    count[k]=len(i)
    
    
plt.figure()
plt.scatter(s,count)
plt.plot(diff(s))
plt.show()