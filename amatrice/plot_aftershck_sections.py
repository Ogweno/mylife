from numpy import genfromtxt,array,zeros,linspace,where,r_,unique,tile
from matplotlib import pyplot as plt
from pyproj import Geod

#afters=genfromtxt('/Users/dmelgar/Amatrice2016/afters/afters.xyz')
#afters_x=afters[:,0]
#afters_y=afters[:,1]


afters_load=genfromtxt('/Users/dmelgar/Amatrice2016/afters/nonlinloc.dat')
afters=zeros((len(afters_load),3))
afters[:,0]=afters_load[:,16]
afters[:,1]=afters_load[:,14]
afters[:,2]=afters_load[:,18]
afters_x=afters_load[:,16]
afters_y=afters_load[:,14]

dist_profile=0.035 #in degs

#Profiles
AAp=array([[12.7575094963,42.9344070448],[13.2606118885,43.082797542]])
BBp=array([[12.8007434665,42.8491552133],[13.2981912131,42.9996964758]])
CCp=array([[12.8449467564,42.7623767325],[13.3402817039,42.9125454083]])
DDp=array([[12.8924346965,42.6819127977],[13.3881579543,42.8261400386]])
EEp=array([[12.9378147046,42.5963626755],[13.4330185746,42.7454790381]])
FFp=array([[12.9859751788,42.5114550687],[13.4729123304,42.6626175681]])

#Make lines
Npts=100
AA=zeros((Npts,2))
BB=zeros((Npts,2))
CC=zeros((Npts,2))
DD=zeros((Npts,2))
EE=zeros((Npts,2))
FF=zeros((Npts,2))
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

#find events within profile "influence"
for kpoint in range(Npts):
    dist=((afters_x-AA[kpoint,0])**2+(afters_y-AA[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
AAafters=afters[keep,0:3]

for kpoint in range(Npts):
    dist=((afters_x-BB[kpoint,0])**2+(afters_y-BB[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
BBafters=afters[keep,0:3]

for kpoint in range(Npts):
    dist=((afters_x-CC[kpoint,0])**2+(afters_y-CC[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
CCafters=afters[keep,0:3]

for kpoint in range(Npts):
    dist=((afters_x-DD[kpoint,0])**2+(afters_y-DD[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
DDafters=afters[keep,0:3]

for kpoint in range(Npts):
    dist=((afters_x-EE[kpoint,0])**2+(afters_y-EE[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
EEafters=afters[keep,0:3]

for kpoint in range(Npts):
    dist=((afters_x-FF[kpoint,0])**2+(afters_y-FF[kpoint,1])**2)**0.5
    i=where(dist<dist_profile)[0]
    if kpoint==0:
        keep=i.copy()
    else:
        keep=r_[keep,i]
keep=unique(keep)
FFafters=afters[keep,0:3]
 
#Get distance from edge of profile
p=Geod(ellps='WGS84')

ref_point=tile(AA[0,:],(len(AAafters),1))
az,baz,dist=p.inv(AAafters[:,0],AAafters[:,1],ref_point[:,0],ref_point[:,1])
AAkm=dist/1000.

ref_point=tile(BB[0,:],(len(BBafters),1))
az,baz,dist=p.inv(BBafters[:,0],BBafters[:,1],ref_point[:,0],ref_point[:,1])
BBkm=dist/1000.

ref_point=tile(CC[0,:],(len(CCafters),1))
az,baz,dist=p.inv(CCafters[:,0],CCafters[:,1],ref_point[:,0],ref_point[:,1])
CCkm=dist/1000.

ref_point=tile(DD[0,:],(len(DDafters),1))
az,baz,dist=p.inv(DDafters[:,0],DDafters[:,1],ref_point[:,0],ref_point[:,1])
DDkm=dist/1000.

ref_point=tile(EE[0,:],(len(EEafters),1))
az,baz,dist=p.inv(EEafters[:,0],EEafters[:,1],ref_point[:,0],ref_point[:,1])
EEkm=dist/1000.

ref_point=tile(FF[0,:],(len(FFafters),1))
az,baz,dist=p.inv(FFafters[:,0],FFafters[:,1],ref_point[:,0],ref_point[:,1])
FFkm=dist/1000.
    
    
#MAKE PLOTS
plt.figure()
plt.scatter(afters[:,0],afters[:,1],s=20,lw=0,c='#808080')
plt.scatter(AAafters[:,0],AAafters[:,1],c='r')
plt.scatter(BBafters[:,0],BBafters[:,1],c='g')
plt.scatter(CCafters[:,0],CCafters[:,1],c='b')
plt.scatter(DDafters[:,0],DDafters[:,1],c='m')
plt.scatter(EEafters[:,0],EEafters[:,1],c='k')
plt.scatter(FFafters[:,0],FFafters[:,1],c='y')
plt.plot(AA[:,0],AA[:,1],'k')
plt.plot(BB[:,0],BB[:,1],'k')
plt.plot(CC[:,0],CC[:,1],'k')
plt.plot(DD[:,0],DD[:,1],'k')
plt.plot(EE[:,0],EE[:,1],'k')
plt.plot(FF[:,0],FF[:,1],'k')
plt.axis('equal')

plt.figure(figsize=(4,20))
plt.subplot(611)
plt.scatter(AAkm,-AAafters[:,2],c='r')
plt.ylabel('Depth (km)')
plt.axis('equal')
plt.subplot(612)
plt.scatter(BBkm,-BBafters[:,2],c='g')
plt.ylabel('Depth (km)')
plt.axis('equal')
plt.subplot(613)
plt.scatter(CCkm,-CCafters[:,2],c='b')
plt.ylabel('Depth (km)')
plt.axis('equal')
plt.subplot(614)
plt.scatter(DDkm,-DDafters[:,2],c='m')
plt.ylabel('Depth (km)')
plt.axis('equal')
plt.subplot(615)
plt.scatter(EEkm,-EEafters[:,2],c='k')
plt.ylabel('Depth (km)')
plt.axis('equal')
plt.subplot(616)
plt.scatter(FFkm,-FFafters[:,2],c='y')
plt.ylabel('Depth (km)')
plt.xlabel('Distance along profile (km)')
plt.axis('equal')
plt.subplots_adjust(left=0.19)



plt.show()
