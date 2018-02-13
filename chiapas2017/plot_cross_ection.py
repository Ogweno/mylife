from matplotlib import pyplot as plt
from numpy import genfromtxt,where,isnan,argmax,ones,array,r_
from pyproj import Geod

bathy=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA_bathy.txt')
slab=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA_slab.txt')
afters=genfromtxt('/Users/dmelgar/Chiapas2017/afters/afters_6days_clean.txt')
AA=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA.xy')
afters_dist=30.0*1000

#Find trench
i=where(isnan(slab[:,2])==False)[0]
slab=slab[i,:]
slab[:,0]=slab[:,0]-360

i=argmax(slab[:,2])
trench=slab[i,0:2]

#Get distances from bathy  to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(bathy[:,0],bathy[:,1],trench[0]*ones(len(bathy)),trench[1]*ones(len(bathy)))
dist=dist/1000
i=where(bathy[:,1]<trench[1])
dist[i]=-dist[i]
dist_bathy=dist
bathy=bathy[:,2]/1000


#Get distances from slab  to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(slab[:,0],slab[:,1],trench[0]*ones(len(slab)),trench[1]*ones(len(slab)))
dist=dist/1000
i=where(slab[:,1]<trench[1])
dist[i]=-dist[i]
dist_slab=dist
slab=slab[:,2]

#Get aftershocks within some distance of profile
af=array([])
count=0
for k in range(len(afters)):
    az,baz,dist=p.inv(AA[:,0],AA[:,1],afters[k,4]*ones(len(AA)),afters[k,3]*ones(len(AA)))
    if dist.min()<afters_dist:
        if count==0:
            af=afters[k,2:6]
            count+=1
        else:    
            af=r_[af,afters[k,2:6]]
af=af.reshape(len(af)/4,4)
af[:,3]=-af[:,3]

#Distance to trench
az,baz,dist=p.inv(af[:,2],af[:,1],trench[0]*ones(len(af)),trench[1]*ones(len(af)))
dist=dist/1000
i=where(af[:,1]<trench[1])
dist[i]=-dist[i]
dist_af=dist


#SSN hypo to trench
#-94.11,14.85,58
az,baz,dist=p.inv(-94.11,14.85,trench[0],trench[1])
dist=dist/1000
dist_ssn=dist



plt.figure()
plt.scatter(dist_af,af[:,3],c='r',s=20)
plt.scatter(dist_ssn,-58,marker='*',s=60)
plt.legend(['Aftershocks','SSN Hypo'])

plt.plot(dist_bathy,bathy,'k')
plt.plot(dist_slab,slab,'k')


plt.xlabel('Distance from trench (km)')
plt.ylabel('Depth (km)')


plt.show()

