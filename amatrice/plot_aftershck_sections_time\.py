from numpy import genfromtxt,array,zeros,linspace,where,r_,unique,tile,fliplr,savetxt,c_
from matplotlib import pyplot as plt
from pyproj import Geod
from scipy.linalg import norm
from obspy.core import UTCDateTime
from string import rjust

min_mag=0.5
dist_profile=0.023 #in degs
dist_profile_faults=0.015 #in degs
yl=[-15,0]
xl=[10,40]
msize=9
Nslices=30


path=u'/Users/dmelgar/code/GMT/Amatrice/'
af1=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af2=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records2.txt',usecols=[0,1,2,3,4,5])
af1_date=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1],dtype='S')
af2_date=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records2.txt',usecols=[0,1],dtype='S')

vettore=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_main.fault')
vettore_synth=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_complex/vettore_synth.fault')
pian=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_complex/pian_grande.fault')
anti=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_listric/antithetic_listric.fault')
detach=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_listric/detachment_large.fault')

def projection_distances(ref_line,aftershock_locations):
    p=Geod(ellps='WGS84')
    dist=zeros(len(aftershock_locations))
    for k in range(len(aftershock_locations)):
        unit_vector=ref_line[-1,:]-ref_line[0,:]
        unit_vector=unit_vector/norm(unit_vector)
        point=aftershock_locations[k,:]
        point=point-ref_line[0,:]
        projected_point=unit_vector.dot(point)*unit_vector+ref_line[0,:]
        az,baz,dist[k]=p.inv(projected_point[0],projected_point[1],ref_line[0,0],ref_line[0,1])
    return dist


i=where(af1[:,5]>min_mag)[0]
j=where(af2[:,5]>min_mag)[0]
print 'Selceted %d events' % (len(i))

afters=r_[af1[i,:],af2[j,:]]

afters_x=af1[i,3]
afters_y=af1[i,2]
afters_t=af1_date
zs1=af1[i,4]
afters_x=r_[afters_x,af2[j,3]]
afters_y=r_[afters_y,af2[j,2]]
afters_t=r_[afters_t,af2_date]
zs2=af2[j,4]

#Convert times to days since origin
for k in range(len(afters_t)):
    time=afters_t[k,0]+'T'+afters_t[k,1]
    time=time.replace('/','-')
    if k==0:
        t0=UTCDateTime(time)
        time_out=array([UTCDateTime(time)-t0])
    else:
        time_out=r_[time_out,array([UTCDateTime(time)-t0])]
afters_t=time_out

#Profiles
profiles=array([[ 13.01567151,  42.54648126,  13.51355551,  42.71660122],
       [ 12.98991192,  42.58726628,  13.48812229,  42.75738348],
       [ 12.96411872,  42.62804522,  13.4626561 ,  42.79815965],
       [ 12.93829183,  42.66881804,  13.43715684,  42.83892972],
       [ 12.91243116,  42.70958475,  13.41162445,  42.87969366],
       [ 12.88653663,  42.75034532,  13.38605883,  42.92045146],
       [ 12.86060814,  42.79109973,  13.36045989,  42.96120311],
       [ 12.83464561,  42.83184798,  13.33482755,  43.00194859],
       [ 12.80864894,  42.87259004,  13.30916172,  43.04268788],
       [ 12.78261806,  42.9133259 ,  13.28346231,  43.08342096]])
AAp=array([profiles[0,[0,2]],profiles[0,[1,3]]]).T
BBp=array([profiles[1,[0,2]],profiles[1,[1,3]]]).T
CCp=array([profiles[2,[0,2]],profiles[2,[1,3]]]).T
DDp=array([profiles[3,[0,2]],profiles[3,[1,3]]]).T
EEp=array([profiles[4,[0,2]],profiles[4,[1,3]]]).T
FFp=array([profiles[5,[0,2]],profiles[5,[1,3]]]).T
GGp=array([profiles[6,[0,2]],profiles[6,[1,3]]]).T
HHp=array([profiles[7,[0,2]],profiles[7,[1,3]]]).T
IIp=array([profiles[8,[0,2]],profiles[8,[1,3]]]).T
JJp=array([profiles[9,[0,2]],profiles[9,[1,3]]]).T

#Make lines
Npts=100
AA=zeros((Npts,2))
BB=zeros((Npts,2))
CC=zeros((Npts,2))
DD=zeros((Npts,2))
EE=zeros((Npts,2))
FF=zeros((Npts,2))
GG=zeros((Npts,2))
HH=zeros((Npts,2))
II=zeros((Npts,2))
JJ=zeros((Npts,2))
AA[:,0]=linspace(AAp[0,0],AAp[1,0],Npts)
AA[:,1]=linspace(AAp[0,1],AAp[1,1],Npts)
BB[:,0]=linspace(BBp[0,0],BBp[1,0],Npts)
BB[:,1]=linspace(BBp[0,1],BBp[1,1],Npts)
CC[:,0]=linspace(CCp[0,0],CCp[1,0],Npts)
CC[:,1]=linspace(CCp[0,1],CCp[1,1],Npts)
DD[:,0]=linspace(DDp[0,0],DDp[1,0],Npts)
DD[:,1]=linspace(DDp[0,1],DDp[1,1],Npts)
EE[:,0]=linspace(EEp[0,0],EEp[1,0],Npts)
EE[:,1]=linspace(EEp[0,1],EEp[1,1],Npts)
FF[:,0]=linspace(FFp[0,0],FFp[1,0],Npts)
FF[:,1]=linspace(FFp[0,1],FFp[1,1],Npts)
GG[:,0]=linspace(GGp[0,0],GGp[1,0],Npts)
GG[:,1]=linspace(GGp[0,1],GGp[1,1],Npts)
HH[:,0]=linspace(HHp[0,0],HHp[1,0],Npts)
HH[:,1]=linspace(HHp[0,1],HHp[1,1],Npts)
II[:,0]=linspace(IIp[0,0],IIp[1,0],Npts)
II[:,1]=linspace(IIp[0,1],IIp[1,1],Npts)
JJ[:,0]=linspace(JJp[0,0],JJp[1,0],Npts)
JJ[:,1]=linspace(JJp[0,1],JJp[1,1],Npts)

#Find vettore faults within profile
def get_fault_points(fault,profile,dist_profile_faults,Npts):
    for kpoint in range(Npts):
        dist=((fault[:,1]-profile[kpoint,0])**2+(fault[:,2]-profile[kpoint,1])**2)**0.5
        i=where(dist<dist_profile_faults)[0]
        if kpoint==0:
            keep=i.copy()
        else:
            keep=r_[keep,i]
    keep=unique(keep)
    points=fault[keep,1:4]
    dist=projection_distances(profile,points[:,0:2])
    km_dist=dist/1000.
    return km_dist,points

AAkm_vettore,AAvettore=get_fault_points(vettore,AA,dist_profile_faults,Npts)
BBkm_vettore,BBvettore=get_fault_points(vettore,BB,dist_profile_faults,Npts)
CCkm_vettore,CCvettore=get_fault_points(vettore,CC,dist_profile_faults,Npts)
DDkm_vettore,DDvettore=get_fault_points(vettore,DD,dist_profile_faults,Npts)
EEkm_vettore,EEvettore=get_fault_points(vettore,EE,dist_profile_faults,Npts)
FFkm_vettore,FFvettore=get_fault_points(vettore,FF,dist_profile_faults,Npts)
GGkm_vettore,GGvettore=get_fault_points(vettore,GG,dist_profile_faults,Npts)
HHkm_vettore,HHvettore=get_fault_points(vettore,HH,dist_profile_faults,Npts)
IIkm_vettore,IIvettore=get_fault_points(vettore,II,dist_profile_faults,Npts)
JJkm_vettore,JJvettore=get_fault_points(vettore,JJ,dist_profile_faults,Npts)

AAkm_vettore_synth,AAvettore_synth=get_fault_points(vettore_synth,AA,dist_profile_faults,Npts)
BBkm_vettore_synth,BBvettore_synth=get_fault_points(vettore_synth,BB,dist_profile_faults,Npts)
CCkm_vettore_synth,CCvettore_synth=get_fault_points(vettore_synth,CC,dist_profile_faults,Npts)
DDkm_vettore_synth,DDvettore_synth=get_fault_points(vettore_synth,DD,dist_profile_faults,Npts)
EEkm_vettore_synth,EEvettore_synth=get_fault_points(vettore_synth,EE,dist_profile_faults,Npts)
FFkm_vettore_synth,FFvettore_synth=get_fault_points(vettore_synth,FF,dist_profile_faults,Npts)
GGkm_vettore_synth,GGvettore_synth=get_fault_points(vettore_synth,GG,dist_profile_faults,Npts)
HHkm_vettore_synth,HHvettore_synth=get_fault_points(vettore_synth,HH,dist_profile_faults,Npts)
IIkm_vettore_synth,IIvettore_synth=get_fault_points(vettore_synth,II,dist_profile_faults,Npts)
JJkm_vettore_synth,JJvettore_synth=get_fault_points(vettore_synth,JJ,dist_profile_faults,Npts)

AAkm_pian,AApian=get_fault_points(pian,AA,dist_profile_faults,Npts)
BBkm_pian,BBpian=get_fault_points(pian,BB,dist_profile_faults,Npts)
CCkm_pian,CCpian=get_fault_points(pian,CC,dist_profile_faults,Npts)
DDkm_pian,DDpian=get_fault_points(pian,DD,dist_profile_faults,Npts)
EEkm_pian,EEpian=get_fault_points(pian,EE,dist_profile_faults,Npts)
FFkm_pian,FFpian=get_fault_points(pian,FF,dist_profile_faults,Npts)
GGkm_pian,GGpian=get_fault_points(pian,GG,dist_profile_faults,Npts)
HHkm_pian,HHpian=get_fault_points(pian,HH,dist_profile_faults,Npts)
IIkm_pian,IIpian=get_fault_points(pian,II,dist_profile_faults,Npts)
JJkm_pian,JJpian=get_fault_points(pian,JJ,dist_profile_faults,Npts)

AAkm_anti,AAanti=get_fault_points(anti,AA,dist_profile_faults,Npts)
BBkm_anti,BBanti=get_fault_points(anti,BB,dist_profile_faults,Npts)
CCkm_anti,CCanti=get_fault_points(anti,CC,dist_profile_faults,Npts)
DDkm_anti,DDanti=get_fault_points(anti,DD,dist_profile_faults,Npts)
EEkm_anti,EEanti=get_fault_points(anti,EE,dist_profile_faults,Npts)
FFkm_anti,FFanti=get_fault_points(anti,FF,dist_profile_faults,Npts)
GGkm_anti,GGanti=get_fault_points(anti,GG,dist_profile_faults,Npts)
HHkm_anti,HHanti=get_fault_points(anti,HH,dist_profile_faults,Npts)
IIkm_anti,IIanti=get_fault_points(anti,II,dist_profile_faults,Npts)
JJkm_anti,JJanti=get_fault_points(anti,JJ,dist_profile_faults,Npts)

AAkm_detach,AAdetach=get_fault_points(detach,AA,dist_profile_faults,Npts)
BBkm_detach,BBdetach=get_fault_points(detach,BB,dist_profile_faults,Npts)
CCkm_detach,CCdetach=get_fault_points(detach,CC,dist_profile_faults,Npts)
DDkm_detach,DDdetach=get_fault_points(detach,DD,dist_profile_faults,Npts)
EEkm_detach,EEdetach=get_fault_points(detach,EE,dist_profile_faults,Npts)
FFkm_detach,FFdetach=get_fault_points(detach,FF,dist_profile_faults,Npts)
GGkm_detach,GGdetach=get_fault_points(detach,GG,dist_profile_faults,Npts)
HHkm_detach,HHdetach=get_fault_points(detach,HH,dist_profile_faults,Npts)
IIkm_detach,IIdetach=get_fault_points(detach,II,dist_profile_faults,Npts)
JJkm_detach,JJdetach=get_fault_points(detach,JJ,dist_profile_faults,Npts)





#find events within profile "influence"
for kpoint in range(Npts):
    dist=((afters_x-AA[kpoint,0])**2+(afters_y-AA[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
AAafters=afters[keep,2:5]
AAafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-BB[kpoint,0])**2+(afters_y-BB[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
BBafters=afters[keep,2:5]
BBafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-CC[kpoint,0])**2+(afters_y-CC[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
CCafters=afters[keep,2:5]
CCafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-DD[kpoint,0])**2+(afters_y-DD[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
DDafters=afters[keep,2:5]
DDafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-EE[kpoint,0])**2+(afters_y-EE[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
EEafters=afters[keep,2:5]
EEafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-FF[kpoint,0])**2+(afters_y-FF[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
FFafters=afters[keep,2:5]
FFafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-GG[kpoint,0])**2+(afters_y-GG[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
GGafters=afters[keep,2:5]
GGafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-HH[kpoint,0])**2+(afters_y-HH[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
HHafters=afters[keep,2:5]
HHafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-II[kpoint,0])**2+(afters_y-II[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
IIafters=afters[keep,2:5]
IIafters_time=afters_t[keep]

for kpoint in range(Npts):
    dist=((afters_x-JJ[kpoint,0])**2+(afters_y-JJ[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
JJafters=afters[keep,2:5]
JJafters_time=afters_t[keep]

#Save events
savetxt(path+'AA.txt',c_[AAafters[:,1],AAafters[:,0]],fmt='%.4f')
savetxt(path+'BB.txt',c_[BBafters[:,1],BBafters[:,0]],fmt='%.4f')
savetxt(path+'CC.txt',c_[CCafters[:,1],CCafters[:,0]],fmt='%.4f')
savetxt(path+'DD.txt',c_[DDafters[:,1],DDafters[:,0]],fmt='%.4f')
savetxt(path+'EE.txt',c_[EEafters[:,1],EEafters[:,0]],fmt='%.4f')
savetxt(path+'FF.txt',c_[FFafters[:,1],FFafters[:,0]],fmt='%.4f')
savetxt(path+'GG.txt',c_[GGafters[:,1],GGafters[:,0]],fmt='%.4f')
savetxt(path+'HH.txt',c_[HHafters[:,1],HHafters[:,0]],fmt='%.4f')
savetxt(path+'II.txt',c_[IIafters[:,1],IIafters[:,0]],fmt='%.4f')
savetxt(path+'JJ.txt',c_[JJafters[:,1],JJafters[:,0]],fmt='%.4f')
 
dist=projection_distances(AA,fliplr(AAafters[:,0:2]))
AAkm=dist/1000.

dist=projection_distances(BB,fliplr(BBafters[:,0:2]))
BBkm=dist/1000.

dist=projection_distances(CC,fliplr(CCafters[:,0:2]))
CCkm=dist/1000.

dist=projection_distances(DD,fliplr(DDafters[:,0:2]))
DDkm=dist/1000.

dist=projection_distances(EE,fliplr(EEafters[:,0:2]))
EEkm=dist/1000.

dist=projection_distances(FF,fliplr(FFafters[:,0:2]))
FFkm=dist/1000.

dist=projection_distances(GG,fliplr(GGafters[:,0:2]))
GGkm=dist/1000.

dist=projection_distances(HH,fliplr(HHafters[:,0:2]))
HHkm=dist/1000.

dist=projection_distances(II,fliplr(IIafters[:,0:2]))
IIkm=dist/1000.

dist=projection_distances(JJ,fliplr(JJafters[:,0:2]))
JJkm=dist/1000.
 
    
#Save positions
savetxt('/Users/dmelgar/Amatrice2016/plots/CC.txt',c_[CCkm,-CCafters[:,2 ]],fmt='%.4f')  
savetxt('/Users/dmelgar/Amatrice2016/plots/CCvettore.txt',c_[CCkm_vettore,-CCvettore[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/CCvettore_synth.txt',c_[CCkm_vettore_synth,-CCvettore_synth[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/CCpian.txt',c_[CCkm_pian,-CCpian[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/CCdetach.txt',c_[CCkm_detach,-CCdetach[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/EE.txt',c_[EEkm,-EEafters[:,2 ]],fmt='%.4f')  
savetxt('/Users/dmelgar/Amatrice2016/plots/EEvettore.txt',c_[EEkm_vettore,-EEvettore[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/EEvettore_synth.txt',c_[EEkm_vettore_synth,-EEvettore_synth[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/EEpian.txt',c_[EEkm_pian,-EEpian[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/EEdetach.txt',c_[EEkm_detach,-EEdetach[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/JJ.txt',c_[JJkm,-JJafters[:,2 ]],fmt='%.4f')  
savetxt('/Users/dmelgar/Amatrice2016/plots/JJvettore.txt',c_[JJkm_vettore,-JJvettore[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/JJvettore_synth.txt',c_[JJkm_vettore_synth,-JJvettore_synth[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/JJpian.txt',c_[JJkm_pian,-JJpian[:,2 ]],fmt='%.4f')
savetxt('/Users/dmelgar/Amatrice2016/plots/JJdetach.txt',c_[JJkm_detach,-JJdetach[:,2 ]],fmt='%.4f')
    
#MAKE PLOTS
dt=afters_t.max()/Nslices
for k in range(Nslices):
    plt.figure(figsize=(8,13))
    t1=k*dt
    t2=(k+1)*dt
    
    i=where((AAafters_time > k*dt) & (AAafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,10)
    if len(i)>0:
        plt.scatter(AAkm[i],-AAafters[i,2],s=msize,lw=0,c='r')
        plt.scatter(AAkm_vettore,-AAvettore[:,2],marker='x',c='k')
        plt.scatter(AAkm_vettore_synth,-AAvettore_synth[:,2],marker='^',c='k')
        plt.scatter(AAkm_anti,-AAanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(AAkm_pian,-AApian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(AAkm_detach,-AAdetach[:,2],marker='s',c='w',s=25,lw=1)
    plt.xlabel('Distance along profile (km)')
    ax.set_yticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    i=where((BBafters_time > k*dt) & (BBafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,9)
    if len(i)>0:
        plt.scatter(BBkm[i],-BBafters[i,2],s=msize,lw=0,c='g')
        plt.scatter(BBkm_vettore,-BBvettore[:,2],marker='x',c='k')
        plt.scatter(BBkm_vettore_synth,-BBvettore_synth[:,2],marker='^',c='k')
        plt.scatter(BBkm_anti,-BBanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(BBkm_pian,-BBpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(BBkm_detach,-BBdetach[:,2],marker='s',c='w',s=25,lw=1)
    plt.ylabel('Depth (km)')
    plt.xlabel('Distance along profile (km)')
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((CCafters_time > k*dt) & (CCafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,8)
    if len(i)>0:
        plt.scatter(CCkm[i],-CCafters[i,2],s=msize,lw=0,c='b')
        plt.scatter(CCkm_vettore,-CCvettore[:,2],marker='x',c='k')
        plt.scatter(CCkm_vettore_synth,-CCvettore_synth[:,2],marker='^',c='k')
        plt.scatter(CCkm_anti,-CCanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(CCkm_pian,-CCpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(CCkm_detach,-CCdetach[:,2],marker='s',c='w',s=25,lw=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((DDafters_time > k*dt) & (DDafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,7)
    if len(i)>0:
        plt.scatter(DDkm[i],-DDafters[i,2],s=msize,lw=0,c='m')
        plt.scatter(DDkm_vettore,-DDvettore[:,2],marker='x',c='k')
        plt.scatter(DDkm_vettore_synth,-DDvettore_synth[:,2],marker='^',c='k')
        plt.scatter(DDkm_anti,-DDanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(DDkm_pian,-DDpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(DDkm_detach,-DDdetach[:,2],marker='s',c='w',s=25,lw=1)
    plt.ylabel('Depth (km)')
    ax.set_xticklabels([])   
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((EEafters_time > k*dt) & (EEafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,6)
    if len(i)>0:
        plt.scatter(EEkm[i],-EEafters[i,2],s=msize,lw=0,c='k')
        plt.scatter(EEkm_vettore,-EEvettore[:,2],marker='x',c='r')
        plt.scatter(EEkm_vettore_synth,-EEvettore_synth[:,2],marker='^',c='r')
        plt.scatter(EEkm_anti,-EEanti[:,2],marker='o',c='R',s=20,lw=1)
        plt.scatter(EEkm_pian,-EEpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(EEkm_detach,-EEdetach[:,2],marker='s',c='w',s=25,lw=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((FFafters_time > k*dt) & (FFafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,5)
    if len(i)>0:
        plt.scatter(FFkm[i],-FFafters[i,2],s=msize,lw=0,c='y')
        plt.scatter(FFkm_vettore,-FFvettore[:,2],marker='x',c='k')
        plt.scatter(FFkm_vettore_synth,-FFvettore_synth[:,2],marker='^',c='k')
        plt.scatter(FFkm_anti,-FFanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(FFkm_pian,-FFpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(FFkm_detach,-FFdetach[:,2],marker='s',c='w',s=25,lw=1)
    plt.ylabel('Depth (km)')
    ax.set_xticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((GGafters_time > k*dt) & (GGafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,4)
    if len(i)>0:
        plt.scatter(GGkm[i],-GGafters[i,2],s=msize,lw=0,c='r')
        plt.scatter(GGkm_vettore,-GGvettore[:,2],marker='x',c='k')
        plt.scatter(GGkm_vettore_synth,-GGvettore_synth[:,2],marker='^',c='k')
        plt.scatter(GGkm_anti,-GGanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(GGkm_pian,-GGpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(GGkm_detach,-GGdetach[:,2],marker='s',c='w',s=25,lw=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((HHafters_time > k*dt) & (HHafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,3)
    if len(i)>0:
        plt.scatter(HHkm[i],-HHafters[i,2],s=msize,lw=0,c='g')
        plt.scatter(HHkm_vettore,-HHvettore[:,2],marker='x',c='k')
        plt.scatter(HHkm_vettore_synth,-HHvettore_synth[:,2],marker='^',c='k')
        plt.scatter(HHkm_anti,-HHanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(HHkm_pian,-HHpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(HHkm_detach,-HHdetach[:,2],marker='s',c='w',s=25,lw=1)
        plt.ylabel('Depth (km)')
    ax.set_xticklabels([])    
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((IIafters_time > k*dt) & (IIafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,2)
    if len(i)>0:
        plt.scatter(IIkm[i],-IIafters[i,2],s=msize,lw=0,c='b')
        plt.scatter(IIkm_vettore,-IIvettore[:,2],marker='x',c='k')
        plt.scatter(IIkm_vettore_synth,-IIvettore_synth[:,2],marker='^',c='k')
        plt.scatter(IIkm_anti,-IIanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(IIkm_pian,-IIpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(IIkm_detach,-IIdetach[:,2],marker='s',c='w',s=25,lw=1)
    ax.set_yticklabels([])
    ax.set_xticklabels([]) 
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    i=where((JJafters_time > k*dt) & (JJafters_time <= (k+1)*dt))[0]
    ax=plt.subplot(5,2,1)
    if len(i)>0:
        plt.scatter(JJkm[i],-JJafters[i,2],s=msize,lw=0,c='m')
        plt.scatter(JJkm_vettore,-JJvettore[:,2],marker='x',c='k')
        plt.scatter(JJkm_vettore_synth,-JJvettore_synth[:,2],marker='^',c='k')
        plt.scatter(JJkm_anti,-JJanti[:,2],marker='o',c='w',s=20,lw=1)
        plt.scatter(JJkm_pian,-JJpian[:,2],marker='>',c='w',s=20,lw=1)
        plt.scatter(JJkm_detach,-JJdetach[:,2],marker='s',c='w',s=25,lw=1)
    plt.ylabel('Depth (km)')
    ax.set_xticklabels([])
    plt.ylim(yl)
    plt.xlim(xl)
    
    
    plt.subplots_adjust(left=0.1,bottom=0.05,top=0.96,right=0.9,hspace=0.05,wspace=0.02)
    plt.suptitle('Events between day %.2f and day %.2f after 8/24' %(t1/86400.,t2/86400))    
    
    plt.savefig('/Users/dmelgar/Amatrice2016/plots/afters/time/cross_sections_'+rjust(str(k),3,'0')+'.png')
    plt.close()