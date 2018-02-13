from numpy import genfromtxt,where,arange,argmin,array,r_,zeros,squeeze,isnan,c_,ones
from numpy.linalg import lstsq
import matplotlib.path as mplPath
from matplotlib import pyplot as plt


decimate=100
overlap_poly=genfromtxt('/Users/dmelgar/Amatrice2016/InSAR/insar_overlap.txt')
desc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/descending/T95_Des_Italy.lltnde')
asc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/ascending/T117_Asc_Italy.lltnde')

#Read Alos
asc_alos=genfromtxt('/Users/dmelgar/Amatrice2016/InSAR/ascending/ALOS2/alos_ascending.txt')
i=where(isnan(asc_alos[:,5])==False)[0]
asc_alos=asc_alos[i,:]

#Read ALOS polygon
alospoly=genfromtxt('/Users/dmelgar/Amatrice2016/InSAR/ascending/ALOS2/ALOS_path.txt')
alos_path = mplPath.Path(alospoly)

#Select area for correction
i=where(alos_path.contains_points(asc_alos[:,0:2])==False)[0]
alos_correct=asc_alos[i,:]

#Fit ramp
lon=alos_correct[:,0]
lat=alos_correct[:,1]
los=alos_correct[:,5]
G=c_[lon**2,lat**2,lon*lat,ones(len(lon))]
los
m,a,b,c=lstsq(G,los)

#Apply correction
correction=m[0]*asc_alos[:,0]**2 + m[1]*asc_alos[:,1]**2 + m[2]*asc_alos[:,0]*asc_alos[:,1] + ones(len(asc_alos))*m[3]
asc_alos[:,5]=asc_alos[:,5]-correction



#Path for overalp are of all interferos
bbPath = mplPath.Path(overlap_poly)

#Decimate descending
i=where(bbPath.contains_points(desc[:,0:2])==True)[0]
j=arange(0,len(i),decimate)
i=i[j]
Ndesc=len(i)
desc=desc[i,:]

#Select points in Sentinel Ascending
i=where(bbPath.contains_points(asc[:,0:2])==True)[0]
asc=asc[i,:]

#Select points in ALOS ascending
i=where(bbPath.contains_points(asc_alos[:,0:2])==True)[0]
asc_alos=asc_alos[i,:]

vmin=-200
vmax=200
plt.figure(figsize=(16,6))
plt.subplot(131)
plt.scatter(asc[:,0],asc[:,1],c=asc[:,6],lw=0,vmin=vmin,vmax=vmax)
plt.axis('equal')
plt.colorbar()
plt.title('Ascending')

plt.subplot(132)
plt.scatter(desc[:,0],desc[:,1],c=desc[:,6],lw=0,vmin=vmin,vmax=vmax)
plt.axis('equal')
plt.colorbar()
plt.title('Descending')

plt.subplot(133)
plt.scatter(asc_alos[:,0],asc_alos[:,1],c=asc_alos[:,5],lw=0,vmin=vmin,vmax=vmax)
plt.axis('equal')
plt.colorbar()
plt.title('ALOS2  Asc.')

plt.show()

#Init
Nlooks=3
G=zeros((Nlooks*Ndesc,3*Ndesc))
d=zeros((Nlooks*Ndesc,1))

print Ndesc
#Find points on ascending, closest to descending
for k in range(Ndesc):
    if k==0:
        print 'Working on pixel No:'
    if k%100==0:
        print '..'+str(k)
    
    lon=desc[k,0]
    lat=desc[k,1]
    
    los_desc=desc[k,6]
    lookE_desc=desc[k,3]
    lookN_desc=desc[k,4]
    lookU_desc=desc[k,5]
    
    #Find closest point on Sentinel ascending
    dist=((asc[:,0]-lon)**2+(asc[:,1]-lat)**2)**0.5
    iasc=argmin(dist)
    
    los_asc=asc[iasc,6]
    lookE_asc=asc[iasc,3]
    lookN_asc=asc[iasc,4]
    lookU_asc=asc[iasc,5]
    
    #Find closest point on ALOS2 ascending
    dist=((asc_alos[:,0]-lon)**2+(asc_alos[:,1]-lat)**2)**0.5
    iasc=argmin(dist)
    
    los_asc_alos=asc_alos[iasc,5]
    lookE_asc_alos=asc_alos[iasc,2]
    lookN_asc_alos=asc_alos[iasc,3]
    lookU_asc_alos=asc_alos[iasc,4]
    
    #Form inversion quantities
    Gtemp=array([[lookE_desc,lookN_desc,lookU_desc],[lookE_asc,lookN_asc,lookU_asc],[lookE_asc_alos,lookN_asc_alos,lookU_asc_alos]])
    dtemp=array([[los_desc],[los_asc],[los_asc_alos]])
    
    G[3*k:3*k+3,3*k:3*k+3]=Gtemp
    d[3*k:3*k+3]=dtemp
        
#Regularization matrix
#L=eye(Ndesc*3)

u,a1,a2,a3=lstsq(G,squeeze(d))
ix=arange(0,len(u),3)
iy=arange(1,len(u),3)
iz=arange(2,len(u),3)
ux=u[ix]      
uy=u[iy]
uz=u[iz]
residual=G.dot(u)-squeeze(d)
res_desc=residual[arange(0,len(residual),3)]
res_asc=residual[arange(1,len(residual),3)]
res_asc_alos=residual[arange(2,len(residual),3)]

plt.figure(figsize=(16,12))
plt.subplot(231)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=ux,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('East (mm)')

plt.subplot(232)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=uy,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('North (mm)')

plt.subplot(233)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=uz,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('Up (mm)')

plt.subplot(234)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=res_desc,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('desc (pred-obs)')

plt.subplot(235)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=res_asc,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('asc (pred-obs)')

plt.subplot(236)
plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=res_asc_alos,lw=0)
plt.axis('equal')
plt.colorbar()
plt.title('asc alos(pred-obs)')