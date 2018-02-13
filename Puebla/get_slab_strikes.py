from numpy import genfromtxt,linspace,where,zeros,savetxt,r_,expand_dims
from scipy.interpolate import interp1d
from pyproj import Geod

slab=genfromtxt('/Users/dmelgar/Slab_Models/ferrari/ferrari_for_strike.txt')

z=[-20,-40,-50,-60,-80]
Npts=50

strike=zeros((Npts*len(z),3))
p=Geod(ellps='WGS84')

#at each depth interpolate and get azimuth between neighboring points
for kz in range(len(z)):
    print kz
    i=where(slab[:,2]==z[kz])[0]
    f=interp1d(slab[i,0],slab[i,1],kind='cubic')
    xout=linspace(slab[i[0],0],slab[i[-1],0],Npts)
    yout=f(xout)
    strike[kz*Npts:(kz+1)*Npts,0]=xout
    strike[kz*Npts:(kz+1)*Npts,1]=yout
    for k in range(kz*Npts,(kz+1)*Npts-1):
        #get azimuth between two points
        az,baz,dist=p.inv(strike[k,0],strike[k,1],strike[k+1,0],strike[k+1,1])
        strike[k,2]=az+180
    #do the final one
    strike[k+1,2]=strike[k,2]
    
    
savetxt(u'/Users/dmelgar/code/GMT/Puebla/strike.xyz',strike,fmt='%.6f')


#make piolygon
i=where(slab[:,2]==-50)[0]
z50=slab[i,0:2]
i=where(slab[:,2]==-80)[0]
z80=slab[i,0:2]
z80=z80[::-1]
out=r_[z50,z80,expand_dims(z50[0,0:2],0)]

savetxt(u'/Users/dmelgar/code/GMT/Puebla/strike_poly.xy',out,fmt='%.6f')