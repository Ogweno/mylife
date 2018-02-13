from numpy import genfromtxt,array,zeros,linspace,where,r_,unique,tile,fliplr,savetxt,c_
from matplotlib import pyplot as plt
from pyproj import Geod
from scipy.linalg import norm

min_mag=0.5
dist_profile=0.025 #profile influence distance in degrees
yl=[-15,0]
xl=[10,40]
msize=5 #marker size used in plots

#load files with aftershocks
path=u'/Users/dmelgar/code/GMT/Amatrice/'
af1=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af2=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records2.txt',usecols=[0,1,2,3,4,5])

#Select aftershcks with a certain minimum magnitude
i=where(af1[:,5]>min_mag)[0]
j=where(af2[:,5]>min_mag)[0]
print 'Selceted %d events' % (len(i))

afters=r_[af1[i,:],af2[j,:]]
afters_x=af1[i,3]
afters_y=af1[i,2]
zs1=af1[i,4]
afters_x=r_[afters_x,af2[j,3]]
afters_y=r_[afters_y,af2[j,2]]
zs2=af2[j,4]


#Define Profile start and end points
profiles=array([[ 12.96411872,  42.62804522,  13.4626561 ,  42.79815965],
       [ 12.91243116,  42.70958475,  13.41162445,  42.87969366],
       [ 12.86060814,  42.79109973,  13.36045989,  42.96120311],
       [ 12.80864894,  42.87259004,  13.30916172,  43.04268788]])

AAp=array([profiles[0,[0,2]],profiles[0,[1,3]]]).T
BBp=array([profiles[1,[0,2]],profiles[1,[1,3]]]).T
CCp=array([profiles[2,[0,2]],profiles[2,[1,3]]]).T
DDp=array([profiles[3,[0,2]],profiles[3,[1,3]]]).T


#Make lines with Npts
Npts=100
AA=zeros((Npts,2))
BB=zeros((Npts,2))
CC=zeros((Npts,2))
DD=zeros((Npts,2))
AA[:,0]=linspace(AAp[0,0],AAp[1,0],Npts)
AA[:,1]=linspace(AAp[0,1],AAp[1,1],Npts)
BB[:,0]=linspace(BBp[0,0],BBp[1,0],Npts)
BB[:,1]=linspace(BBp[0,1],BBp[1,1],Npts)
CC[:,0]=linspace(CCp[0,0],CCp[1,0],Npts)
CC[:,1]=linspace(CCp[0,1],CCp[1,1],Npts)
DD[:,0]=linspace(DDp[0,0],DDp[1,0],Npts)
DD[:,1]=linspace(DDp[0,1],DDp[1,1],Npts)


#This function finds events within profile "influence distance"
def get_afters_profile(afters_x,afters_y,profile_points):
    
    for kpoint in range(Npts):
        dist=((afters_x-profile_points[kpoint,0])**2+(afters_y-profile_points[kpoint,1])**2)**0.5
        i=where(dist<dist_profile)[0]
        if kpoint==0:
            keep=i.copy()
        else:
            keep=r_[keep,i]
    keep=unique(keep)
    selected_aftershocks=afters[keep,2:5]
    return selected_aftershocks

#Select events close toe ach profile
AAafters=get_afters_profile(afters_x,afters_y,AA)
BBafters=get_afters_profile(afters_x,afters_y,BB)
CCafters=get_afters_profile(afters_x,afters_y,CC)
DDafters=get_afters_profile(afters_x,afters_y,DD)

#This function projects the selected aftershocks to "distance along profile"
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

#now convert selected events from coordinates to "distance along profile"
dist=projection_distances(AA,fliplr(AAafters[:,0:2]))
AAkm=dist/1000.

dist=projection_distances(BB,fliplr(BBafters[:,0:2]))
BBkm=dist/1000.

dist=projection_distances(CC,fliplr(CCafters[:,0:2]))
CCkm=dist/1000.

dist=projection_distances(DD,fliplr(DDafters[:,0:2]))
DDkm=dist/1000.


 
    

    
#MAKE PLOTS
plt.figure()
plt.scatter(afters[:,3],afters[:,2],s=20,lw=0,c='#808080')
plt.scatter(AAafters[:,1],AAafters[:,0],c='r')
plt.scatter(BBafters[:,1],BBafters[:,0],c='g')
plt.scatter(CCafters[:,1],CCafters[:,0],c='b')
plt.scatter(DDafters[:,1],DDafters[:,0],c='m')
plt.plot(AA[:,0],AA[:,1],'k')
plt.plot(BB[:,0],BB[:,1],'k')
plt.plot(CC[:,0],CC[:,1],'k')
plt.plot(DD[:,0],DD[:,1],'k')
plt.axis('equal')

plt.figure(figsize=(7,16))
plt.subplot(4,1,4)
plt.scatter(AAkm,-AAafters[:,2],s=msize,lw=0,c='r')
plt.ylabel('Depth (km)')
plt.xlabel('Distance along profile (km)')
plt.ylim(yl)
plt.xlim(xl)

plt.subplot(4,1,3)
plt.scatter(BBkm,-BBafters[:,2],s=msize,lw=0,c='g')
plt.ylabel('Depth (km)')
plt.ylim(yl)
plt.xlim(xl)

plt.subplot(4,1,2)
plt.scatter(CCkm,-CCafters[:,2],s=msize,lw=0,c='b')
plt.ylabel('Depth (km)')
plt.ylim(yl)
plt.xlim(xl)

plt.subplot(4,1,1)
plt.scatter(DDkm,-DDafters[:,2],s=msize,lw=0,c='m')
plt.ylabel('Depth (km)')
plt.ylim(yl)
plt.xlim(xl)

plt.subplots_adjust(left=0.19,bottom=0.05,top=0.99)

plt.show()
