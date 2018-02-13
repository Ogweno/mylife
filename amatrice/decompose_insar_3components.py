from matplotlib import pyplot as plt
from numpy import genfromtxt,argmin,array,zeros,ones,where,linspace,r_,arange,squeeze
from numpy.linalg import lstsq
from matplotlib.ticker import MultipleLocator
import matplotlib.path as mplPath

decimate=150
g=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Oct26-30_combined.txt')
sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Oct26-30_combined.txt',usecols=0,dtype='S')
asc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T117_Italy/T117_Italy.lltnde')
desc=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T22_Italy/T22_Italy.lltnde')
overlap_poly=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/T22_T117_overlap.txt')

#Parse GPS
lon_gps=g[:,1]
lat_gps=g[:,2]
north=g[:,4]
east=g[:,3]
up=g[:,5]

#Path for overalp are of all interferos
bbPath = mplPath.Path(overlap_poly)

#Decimate descending
i=where(bbPath.contains_points(desc[:,0:2])==True)[0]
j=arange(0,len(i),decimate)
i=i[j]
Ndesc=len(i)
desc=desc[i,:]

#Same for ascending
i=where(bbPath.contains_points(asc[:,0:2])==True)[0]
j=arange(0,len(i),decimate)
i=i[j]
Nasc=len(i)
asc=asc[i,:]

vmin=-200
vmax=200
plt.figure(figsize=(12,4.5))
plt.subplot(121)
plt.scatter(asc[:,0],asc[:,1],c=asc[:,6],lw=0,vmin=vmin,vmax=vmax)
plt.colorbar()
plt.title('Ascending')

plt.subplot(122)
plt.scatter(desc[:,0],desc[:,1],c=desc[:,6],lw=0,vmin=vmin,vmax=vmax)
plt.colorbar()
plt.title('Descending')



#Now the inversion
#Init
print "Running inversion"
Nlooks=2
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
    
    #Form inversion quantities
    Gtemp=array([[lookE_desc,lookN_desc,lookU_desc],[lookE_asc,lookN_asc,lookU_asc]])
    dtemp=array([[los_desc],[los_asc]])
    
    G[2*k:2*k+2,3*k:3*k+3]=Gtemp
    d[2*k:2*k+2]=dtemp


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
res_desc=residual[arange(0,len(residual),2)]
res_asc=residual[arange(1,len(residual),2)]


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


plt.figure(figsize=(12,8))
plt.subplot(211)
cb=plt.scatter(desc[:,0],desc[:,1],vmin=vmin,vmax=vmax,c=(ux**2+uy**2)**0.5,lw=0)
plt.quiver(desc[:,0],desc[:,1],ux/(ux**2+uy**2)**0.5,uy/(ux**2+uy**2)**0.5,scale=50)
plt.colorbar(cb)

plt.title('East (mm)')


plt.show()