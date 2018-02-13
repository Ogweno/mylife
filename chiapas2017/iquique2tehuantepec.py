from numpy import genfromtxt
from matplotlib import pyplot as plt
from numpy import cos,sin,array,deg2rad,zeros,argmin,c_,savetxt,arange,meshgrid,where,ones
from mudpy import forward

#read entawai
f=genfromtxt('/Users/dmelgar/Iquique2014/GFZ_model/iquique_gfz.rupt')

#Tehuantepec slab
slab=genfromtxt('/Users/dmelgar/Chiapas2017/slab/chiapas_slab.fault')

dlat2=3.5
dlon2=-1.5

dlat=35
dlon=-25

mult=1.5

#Rotate
angle=deg2rad(-60)
R=array([[cos(angle),-sin(angle)],[sin(angle),cos(angle)]])
xy=zeros((len(f),2))
xy[:,0]=f[:,1]+dlon
xy[:,1]=f[:,2]+dlat
x0=-91
y0=16
xy[:,0]=xy[:,0]-x0
xy[:,1]=xy[:,1]-y0
s=f[:,9]
r=90*ones(len(f))

rxy=xy.dot(R)
xy[:,0]=rxy[:,0]+x0+dlon2
xy[:,1]=rxy[:,1]+y0+dlat2
slab[:,1]=slab[:,1]-360

plt.figure()
plt.subplot(211)
plt.scatter(f[:,3],f[:,2],c=f[:,9],s=90,lw=0,cmap=plt.cm.jet) ; plt.colorbar()

plt.subplot(212)
plt.scatter(xy[:,0],xy[:,1],c=f[:,9],s=90,lw=0,cmap=plt.cm.jet)
plt.scatter(slab[:,1],slab[:,2],s=90,edgecolor='k',facecolor='w',alpha=0.1,lw=1.5)
plt.axis('equal')

plt.show()



#Now asign the slip from the closest subfault
ss=zeros(len(slab))
ds=zeros(len(slab))
for k in range(len(slab)):
    print k
    dist=((slab[k,1]-xy[:,0])**2+(slab[k,2]-xy[:,1])**2)**0.5
    i=argmin(dist)
    dmin=min(dist)
    if dmin<0.1:
        ss[k]=s[i]*cos(deg2rad(r[i]))*mult
        ds[k]=s[i]*sin(deg2rad(r[i]))*mult
    else:
        ss[k]=0
        ds[k]=0
    
plt.figure()
plt.scatter(slab[:,1],slab[:,2],c=(ss**2+ds**2)**0.5,lw=0,s=60,cmap=plt.cm.jet)
plt.colorbar()


#Save it
i=where((ss**2+ds**2)**0.5>0)[0]
out=c_[slab[i,0:8],ss[i],ds[i],slab[i,8:10],zeros((len(i),2))]
fout='/Users/dmelgar/Chiapas2017/tsunami/scenario/chiapas_iquique.rupt'
savetxt(fout,out,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.1f\t%.6f\t%.6f\t%.1f\t%.1f\t%.1f\t%.1f')


plt.show()


#Make dtopo
dtopo='/Users/dmelgar/Chiapas2017/tsunami/scenario/chiapas_iquique.dtopo'
dl=0.02
y=arange(13,17,dl)
x=arange(-97,-91,dl)
#X,Y=meshgrid(x,y)
#x=X.ravel()
#y=Y.ravel()
forward.move_seafloor_okada(fout,dtopo,x,y)
