from numpy import genfromtxt,where,savetxt,c_
from pyproj import Proj

af=genfromtxt('/Users/dmelgar/Iquique2014/Aftershocks/replicas_TA_20140401_20150417',delimiter='|',usecols=[0,1,6,7,8])
slab=genfromtxt('/Users/dmelgar/code/GMT/Iquique/slab.xyz')
slab[:,2]=-slab[:,2]
#Keep only certain time span
i=where((af[:,0]==2014) & ((af[:,1]==4) | (af[:,1]==4)))[0]
afout=c_[af[i,3],af[i,2],af[i,4]]
#Filter by distance to the slab (convert to UTM and calculate distance
p=Proj(proj='utm',ellps='WGS84',zone='19K')
afx,afy=p(afout[:,0],afout[:,1])
slabx,slaby=p(slab[:,0],slab[:,1])
afz=afout[:,2]*-1000
slabz=slab[:,2]*1000
i=[]
distance_threshold=15
for k in range(len(afx)):
    distance=(((afx[k]-slabx)**2+(afy[k]-slaby)**2+(afz[k]-slabz)**2)**0.5).min()/1000
    if distance<distance_threshold:
        i.append(k)
print str(len(i))+'/'+str(len(afout))+' aftershocks within '+str(distance_threshold)+'km of slab'
savetxt('/Users/dmelgar/code/GMT/Iquique/afters.xy',afout[i,:],header=' lon, lat, z(km)',fmt='%.6f\t%.6f\t%.2f')