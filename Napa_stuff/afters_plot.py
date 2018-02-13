from numpy import genfromtxt,zeros,where,meshgrid,arange
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata



def ll2utm(lon,lat):
    import utm
    try:
        x=zeros(lon.shape)
        y=zeros(lon.shape)
        for k in range(len(lon)):  
            u=utm.from_latlon(lat[k],lon[k])
            x[k]=u[0]/1000
            y[k]=u[1]/1000
    except:
        u=utm.from_latlon(lat,lon)
        x=u[0]/1000
        y=u[1]/1000
    return x,y  
        

s1=genfromtxt(u'/Users/dmelgar/code/GMT/Napa/fault_rupture.line1.txt')
s1x=s1[:,0]
s1y=s1[:,1]
s1x,s1y=ll2utm(s1x,s1y)
s2=genfromtxt(u'/Users/dmelgar/code/GMT/Napa/fault_rupture.line2.txt')
s2x=s2[:,0]
s2y=s2[:,1]
s2x,s12=ll2utm(s2x,s2y)
af=genfromtxt(u'/Users/dmelgar/Napa2014/hardebeck_afters.txt',usecols=[1,2,3])
afx=af[:,1]
afy=af[:,0]
afx,afy=ll2utm(afx,afy)
afz=-af[:,2]
#f=genfromtxt('/Users/dmelgar/Slip_inv/Napa_fcmt/data/model_info/fcmt.fault',usecols=[1,2,3])
f=genfromtxt('/Users/dmelgar/Slip_inv/misc/surface_slip.fault',usecols=[1,2,3])
fx=f[:,0]
fy=f[:,1]
fx,fy=ll2utm(fx,fy)
fz=-f[:,2]
f1=arange(fx.min(),fx.max(),0.2)
f2=arange(fy.min(),fy.max(),0.2)
FX,FY=meshgrid(f1,f2)
FZ=griddata(fx,fy,fz,FX,FY)

u=genfromtxt('/Users/dmelgar/Faults/west_napa_ucerf2.txt')
ux=u[:,1]
uy=u[:,0]
ux,uy=ll2utm(ux,uy)

lonmin=-122.4
lonmax=-122.25
latmin=38.0
latmax=38.45
lonmin,latmin=ll2utm(lonmin,latmin)
lonmax,latmax=ll2utm(lonmax,latmax)

#Filter afters
i=where((afx>lonmin) & (afx<lonmax) & (afy>latmin) & (afy<latmax))[0]
afx=afx[i]
afy=afy[i]
afz=afz[i]
#Filter UCERF
i=where((ux>lonmin) & (ux<lonmax) & (uy>latmin) & (uy<latmax))[0]
ux=ux[i]
uy=uy[i]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(s1x,s1y,zeros(s1x.shape),c='k',lw=3)
ax.plot(s2x,s2y,zeros(s2x.shape),c='k',lw=3)
ax.plot(ux,uy,c='brown',lw=3)
ax.scatter(afx,afy,afz)
#ax.scatter(fx,fy,fz,c='g',s=20)
ax.plot_surface(FX,FY,FZ,color='g',shade=False,alpha=0.5,linewidth=0, antialiased=True,rstride=1,cstride=1)
ax.set_xlim(lonmin,lonmax)
ax.set_ylim(latmin,latmax)
ax.set_zlim(-15,0)
ax.set_xlabel('Easting (km)')
ax.set_ylabel('Northing (km)')
ax.set_zlabel('Depth (km)')
plt.show()

