from matplotlib import pyplot as plt
from pyproj import Geod
from numpy import ones,array,genfromtxt,unique,zeros,where,sqrt,c_,savetxt,expand_dims

hypo=array([-98.6878,18.3044,-57.52])
slab=array([[-99.9936,16.3264,-4.485],[-99.6813,16.8054,-20],[-99.3585,17.3107,-40],[-98.5383,18.5437,-60],[-98.4533,18.6852,-80],[-98.3320,18.8620,-100]])
fault=genfromtxt('/Users/dmelgar/Slip_inv/puebla_herloc_NP2/output/inverse_models/models/allsites_II_vr3.2.0016.inv')
topo=genfromtxt(u'/Users/dmelgar/Puebla2017/cross-section/AAp_topo.xyz')

#get moment per depth
z=unique(fault[:,3])
i=where(z>51)
z=z[i]
M0=zeros(len(z))
lon_slip=zeros(len(z))
lat_slip=zeros(len(z))
depth_slip=zeros(len(z))
for k in range(len(z)):
    i=where(fault[:,3]==z[k])[0]
    M0[k]=sum(sqrt(fault[i,8]**2+fault[i,9]**2)*fault[i,10]*fault[i,11]*fault[i,13])
    lon_slip[k]=fault[i[14],1]
    lat_slip[k]=fault[i[14],2]
    depth_slip[k]=fault[i[0],3]



#get distances from trench
p=Geod(ellps='WGS84')
az,baz,slab_dist=p.inv(slab[:,0],slab[:,1],ones(len(slab))*slab[0,0],ones(len(slab))*slab[0,1])
slab_dist=slab_dist/1000

az,baz,hypo_dist=p.inv(hypo[0],hypo[1],slab[0,0],slab[0,1])
hypo_dist=hypo_dist/1000

az,baz,fault_dist=p.inv(lon_slip,lat_slip,ones(len(lon_slip))*slab[0,0],ones(len(lon_slip))*slab[0,1])
fault_dist/=1000

az,baz,topo_dist=p.inv(topo[:,0],topo[:,1],ones(len(topo))*slab[0,0],ones(len(topo))*slab[0,1])
topo_dist/=1000
i=where(topo[:,1]<slab[0,1])[0]
topo_dist[i]=-topo_dist[i]

out=c_[slab_dist+80,slab[:,2]]
savetxt(u'/Users/dmelgar/code/GMT/Puebla/slab.dist',out,fmt='%.4f')

out=c_[topo_dist+80,topo[:,2]/1000]
savetxt(u'/Users/dmelgar/code/GMT/Puebla/topo.dist',out,fmt='%.4f')

out=expand_dims(array([hypo_dist+80,hypo[2]]),0)
savetxt(u'/Users/dmelgar/code/GMT/Puebla/hypo.dist',out,fmt='%.4f')

out=c_[fault_dist+80,-depth_slip,M0/1e18]
savetxt(u'/Users/dmelgar/code/GMT/Puebla/moment.dist',out,fmt='%.4f')
    

plt.figure(figsize=(6,3))
plt.plot(topo_dist,topo[:,2]/1000)
plt.plot(slab_dist,slab[:,2])
plt.plot(slab_dist,slab[:,2]-10)
plt.plot(slab_dist,slab[:,2]-40)
plt.scatter(slab_dist,slab[:,2])
scat=plt.scatter(fault_dist,-depth_slip,marker='o',s=10,c=M0,cmap=plt.cm.jet)
cb=plt.colorbar(scat)
cb.set_label('Moment (Nm)')
plt.scatter(hypo_dist,hypo[2],marker='*',s=80)
plt.grid()
plt.xlim([-70,400])
#plt.axis('equal')

plt.show()