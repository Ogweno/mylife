from matplotlib import pyplot as plt
from numpy import genfromtxt,where,isnan,argmax,ones,array,r_,c_,savetxt,expand_dims,unique,zeros,linspace
from pyproj import Geod
from scipy.interpolate import interp1d

bathy=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA_bathy.txt')
slab=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA_slab.txt')
afters=genfromtxt('/Users/dmelgar/Chiapas2017/afters/SSNMX_cleaned.csv',delimiter=',',usecols=[3,4,5])
AA=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA.xy')
afters_dist=20.0*1000
mts_dist=60.0*1000
thrust_mts_dist=60.0*1000
gmt_path=u'/Users/dmelgar/code/GMT/tehuantepec/'
mts=genfromtxt(u'/Users/dmelgar/Chiapas2017/seismicity/normal_faults.psmeca')
mts_thrust=genfromtxt(u'/Users/dmelgar/Chiapas2017/seismicity/thrust_oaxaca-elsalvador.psmeca')
origin=array([-95.2738,13.4614])
fault=genfromtxt(u'/Users/dmelgar/Slip_inv/Chiapas_hernandez_new/output/inverse_models/models/final_vr_3.4.0007.inv.total')
fault_dist=7.0*1000
slip_thresh=1.0

offset_trench=85.0

#Find trench
i=where(isnan(slab[:,2])==False)[0]
slab=slab[i,:]
slab[:,0]=slab[:,0]-360

i=argmax(slab[:,2])
#trench=slab[i,0:2]
trench=origin
trench2=slab[i,0:2]

#Get distances from bathy  to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(bathy[:,0],bathy[:,1],trench[0]*ones(len(bathy)),trench[1]*ones(len(bathy)))
dist=dist/1000
i=where(bathy[:,1]<trench[1])
dist[i]=-dist[i]
dist_bathy=dist-offset_trench
bathy=bathy[:,2]/1000

#Make seaward portion of bathy
bathy2=genfromtxt('/Users/dmelgar/Chiapas2017/misc/AA_bathy.txt')
az,baz,dist2=p.inv(bathy2[:,0],bathy2[:,1],trench2[0]*ones(len(bathy)),trench2[1]*ones(len(bathy)))
dist2=dist2/1000
i=where(bathy2[:,1]<trench2[1])
dist_bathy_ocean=abs(dist2[i]-dist2[i][0])-offset_trench
bathy_ocean=bathy2[i,2]/1000.


#Get distances from slab  to trench
p=Geod(ellps='WGS84')
az,baz,dist=p.inv(slab[:,0],slab[:,1],trench[0]*ones(len(slab)),trench[1]*ones(len(slab)))
dist=dist/1000
i=where(slab[:,1]<trench[1])
dist[i]=-dist[i]
dist_slab=dist-offset_trench
slab=slab[:,2]


#Get distances  from slab to trench for curvature
curvature=genfromtxt('/Users/dmelgar/Chiapas2017/flexure/AA_slabcurv.txt')
az,baz,dist=p.inv(curvature[:,0],curvature[:,1],trench[0]*ones(len(curvature)),trench[1]*ones(len(curvature)))
dist=dist/1000
i=where(curvature[:,1]<trench[1])
dist[i]=-dist[i]
dist_curvature=dist-offset_trench
curvature=curvature[:,3]*1e6

#Get aftershocks within some distance of profile
af=array([])
count=0
for k in range(len(afters)):
    az,baz,dist=p.inv(AA[:,0],AA[:,1],afters[k,1]*ones(len(AA)),afters[k,0]*ones(len(AA)))
    if dist.min()<afters_dist:
        if count==0:
            af=expand_dims(afters[k,:],0)
            count+=1
        else:    
            af=r_[af,expand_dims(afters[k,:],0)]
af[:,2]=-af[:,2]


#Get mtss within some distance of profile
count=0
for k in range(len(mts)):
    az,baz,dist=p.inv(AA[:,0],AA[:,1],mts[k,0]*ones(len(AA)),mts[k,1]*ones(len(AA)))
    if dist.min()<mts_dist:
        if count==0:
            mt_out=mts[k,:]
            mt_out=expand_dims(mt_out,0)
            #mt_out[0,0]=dist.min()/1000.
            #mt_out[0,1]=mt_out[0,2]
            count+=1
        else:    
            mt=expand_dims(mts[k,:],0)
            #mt[0,0]=dist.min()/1000.
            #mt[0,1]=mt_out[0,2]
            mt_out=r_[mt_out,mt]
mt_out[:,2]=-mt_out[:,2]
savetxt(gmt_path+'normal_mts.psmeca',mt_out,fmt='%.4f')


#Get thrust mtss within some distance of profile
count=0
for k in range(len(mts)):
    az,baz,dist=p.inv(AA[:,0],AA[:,1],mts_thrust[k,0]*ones(len(AA)),mts_thrust[k,1]*ones(len(AA)))
    print dist.min()/1e3
    if dist.min()<thrust_mts_dist:
        if count==0:
            mt_out=mts_thrust[k,:]
            mt_out=expand_dims(mt_out,0)
            #mt_out[0,0]=dist.min()/1000.
            #mt_out[0,1]=mt_out[0,2]
            count+=1
        else:    
            mt=expand_dims(mts_thrust[k,:],0)
            #mt[0,0]=dist.min()/1000.
            #mt[0,1]=mt_out[0,2]
            mt_out=r_[mt_out,mt]
mt_out[:,2]=-mt_out[:,2]
savetxt(gmt_path+'thrust_mts.psmeca',mt_out,fmt='%.4f')


#Get faults within some distance of profile
count=0
for k in range(len(fault)):
    az,baz,dist=p.inv(AA[:,0],AA[:,1],fault[k,1]*ones(len(AA)),fault[k,2]*ones(len(AA)))
    if dist.min()<fault_dist:
        if count==0:
            fault_select=fault[k,1:4]
            fault_select=expand_dims(fault_select,0)
            count+=1
        else:    
            fault_select=r_[fault_select,expand_dims(fault[k,1:4],0)]
fault_select[:,2]=-fault_select[:,2]
fault_out=fault_select[:,1:3]

#Get distance from those faults to profile start
for k in range(len(fault_out)):
    az,baz,dist=p.inv(fault_select[k,0],fault_select[k,1],origin[0],origin[1])
    fault_out[k,0]=dist/1000.-offset_trench



#Now get average slip of faults over a threshold at each depth
unique_depths=unique(fault[:,3])
slip=zeros(len(unique_depths))
for k in range(len(unique_depths)):
    i=where(fault[:,3]==unique_depths[k])[0]
    s=(fault[i,8]**2+fault[i,9]**2)**0.5
    M0=s*10e3*10e3*fault[i,-1]
    slip[k]=M0.sum()

#Interpolate this to something denser
xi=linspace(fault_out[:,0].min(),fault_out[:,0].max(),100)
xi=xi[::-1]
zi=linspace(fault_out[:,1].min(),fault_out[:,1].max(),100)
f=interp1d(fault_out[:,1],slip)
si=f(zi)/1e20
savetxt(gmt_path+'faults_cross.txt',c_[xi,zi,si],fmt='%.4f')


#Distance to trench
az,baz,dist=p.inv(af[:,1],af[:,0],origin[0]*ones(len(af)),origin[1]*ones(len(af)))
dist=dist/1000
dist_af=dist-offset_trench
out=c_[dist_af,af[:,2]]
savetxt(gmt_path+'afters.txt',out,fmt='%.4f')

#SSN hypo to trench
az,baz,dist=p.inv(-94.103,14.761,trench[0],trench[1])
dist=dist/1000-offset_trench
dist_ssn=dist


#Write files that GMT will use
slab=slab+1
out=c_[dist_slab,slab]
savetxt(gmt_path+'AAslab.txt',out,fmt='%.4f')


out=c_[dist_curvature,curvature]
savetxt(gmt_path+'AAcurvature.txt',out,fmt='%.4f')


#Next two are bottom of crust
dist_slab_bot=dist_slab-4
slab_bot=slab-9
out=c_[dist_slab_bot,slab_bot]
savetxt(gmt_path+'AAslab_bot.txt',out,fmt='%.4f')

clip=8
dist_bathy_ocean=dist_bathy_ocean[0:-clip]
bathy_ocean=bathy_ocean[0,0:-clip]-9
out=c_[dist_bathy_ocean,bathy_ocean.T]
savetxt(gmt_path+'AAslab_bot2.txt',out,fmt='%.4f')

#Next 2 are bottom of litho qt 50km depth
clip=11
dist_asth1=dist_bathy_ocean[0:-clip]
asth1=bathy_ocean[0:-clip]-9-29.5
out1=c_[dist_asth1,asth1.T]

dist_asth2=dist_slab-4-7
asth2=slab-9-38
out2=c_[dist_asth2,asth2]
out=r_[out1,out2]
savetxt(gmt_path+'AA_asth.txt',out,fmt='%.4f')

#Next 2 are bottom of litho qt 70km depth
clip=13
dist_asth3=dist_bathy_ocean[0:-clip]
asth3=bathy_ocean[0:-clip]-9-47.5
out1=c_[dist_asth3,asth3.T]

dist_asth4=dist_slab-4-9
asth4=slab-9-55.5
out2=c_[dist_asth4,asth4]
out=r_[out1,out2]
savetxt(gmt_path+'AA_asth70.txt',out,fmt='%.4f')


#Other stuff

out=c_[dist_bathy,bathy]
savetxt(gmt_path+'AAbathy.txt',out,fmt='%.4f')

out=array([[dist_ssn+1],[-46]]).T
savetxt(gmt_path+'AAhypocenter.txt',out,fmt='%.4f')







plt.figure()
plt.scatter(dist_af,af[:,2],c='r',s=20)
plt.scatter(dist_ssn,-46,marker='*',s=60)
plt.legend(['Aftershocks','SSN Hypo'])

plt.plot(dist_bathy,bathy,'k')
plt.plot(dist_slab,slab,'k')
plt.plot(dist_bathy,bathy,'k')
plt.plot(dist_slab_bot,slab_bot,'k')
plt.plot(dist_bathy_ocean,bathy_ocean.squeeze(),'r',lw=1)
plt.plot(dist_asth2,asth2)
plt.plot(dist_asth1,asth1.squeeze(),'r',lw=1)
plt.plot(dist_asth3,asth3)
plt.plot(dist_asth4,asth4.squeeze(),'r',lw=1)
plt.scatter(xi,zi,c=si,cmap=plt.cm.jet)
plt.scatter(fault_out[:,0],fault_out[:,1],c=slip,cmap=plt.cm.jet)
plt.colorbar()
plt.xlim([0,300])

plt.xlabel('Distance from trench (km)')
plt.ylabel('Depth (km)')


print slip.max()
plt.show()

