from numpy import arange,zeros,flipud,fliplr,meshgrid,pi,savetxt,c_,nan,where,r_,c_,ones,isnan,unique,sqrt,argmin
from numpy.linalg import lstsq
import h5py


##Ascending interfero
#f=h5py.File(u'/Users/dmelgar/Lefkada2015/InSAR/S1A_IW_175_3327_20151105-20151117_0011_00014.h5','r')
#out='/Users/dmelgar/Lefkada2015/InSAR/ascending.los'
## f.attrs.items()
#los=flipud(f['GEOCODE'][u'unwrapped_interferogram'].value)
#e=flipud(f['GEOCODE']['line_of_sight_e'])
#n=flipud(f['GEOCODE']['line_of_sight_n'])
#u=flipud(f['GEOCODE']['line_of_sight_u'])
## los.shape (8281, 12601)
#decimate=10
#wavelength=0.055465763000000001/(4*pi)
#delta=0.00027777777777777778
#west=19.699722222222221
#south=37.199722222222221
#east=23.199999999999999
#north=39.5

##Descending interfero
f=h5py.File(u'/Users/dmelgar/Lefkada2015/InSAR/S1A_IW_080_2222_20151111-20151123_0011_00070.h5','r')
out='/Users/dmelgar/Lefkada2015/InSAR/descending.los'
# f.attrs.items()
los=flipud(f['GEOCODE'][u'unwrapped_interferogram'].value)
e=flipud(f['GEOCODE']['line_of_sight_e'])
n=flipud(f['GEOCODE']['line_of_sight_n'])
u=flipud(f['GEOCODE']['line_of_sight_u'])
# los.shape (8281, 12601)
decimate=5
wavelength=0.055465763000000001/(4*pi)
delta=0.00027777777777777778
west=19.300000000000001
south=37.700000000000003
east=22.900000000000002
north=40.0


lon=arange(west,east-delta/2,delta)
lat=arange(south,north-delta,delta)
ilon=arange(0,len(lon)/decimate)*decimate
ilat=arange(0,len(lat)/decimate)*decimate
[X,Y]=meshgrid(lon,lat)


lon_out=zeros(len(ilon)*len(ilat))
lat_out=zeros(len(ilon)*len(ilat))
los_out=zeros(len(ilon)*len(ilat))
lookE=zeros(len(ilon)*len(ilat))
lookN=zeros(len(ilon)*len(ilat))
lookU=zeros(len(ilon)*len(ilat))
k=0
for klon in ilon:
    print str(klon/decimate)+'/'+str(len(ilon))
    for klat in ilat:
        lon_out[k]=X[klat,klon]
        lat_out[k]=Y[klat,klon]
        if los[klat,klon]==0:
            los_out[k]=nan
        else:
            los_out[k]=los[klat,klon]*-wavelength
        lookN[k]=n[klat,klon]
        lookE[k]=e[klat,klon]
        lookU[k]=u[klat,klon]
        k+=1

#Determine points to be used for removing ramp
i=where(lat_out<38.3)[0]
i2=where(lat_out>39)[0]
i3=where((lon_out<20.9) & (lat_out<38.3))[0]
j=where((lon_out>21.15) & (lat_out<39) & (lat_out>38.3))[0]
ij=r_[i,i3,i3,i3,i3,i2,j]
k=where(isnan(los_out[ij])==0)[0]
G=c_[lon_out[ij[k]]**2,lat_out[ij[k]]**2,lon_out[ij[k]]*lat_out[ij[k]],ones(len(ij[k]))]
d=los_out[ij[k]]
m,a,b,c=lstsq(G,d)
los_ou


#Find points clsoes to GPS sites and print los vector
d=sqrt((20.585-lon_out)**2+(38.619-lat_out)**2)
i=argmin(d)
print 'PONT (n,e,u) = ('+str(lookN[i])+','+str(lookE[i])+','+str(lookU[i])+')'
print 'PONT LOS = '+str(los_out[i])
d=sqrt((20.674-lon_out)**2+(38.781-lat_out)**2)
i=argmin(d)
print 'SPAN (n,e,u) = ('+str(lookN[i])+','+str(lookE[i])+','+str(lookU[i])+')'
print 'SPAN LOS = '+str(los_out[i])
d=sqrt((20.589-lon_out)**2+(38.177-lat_out)**2)
i=argmin(d)
print 'VLSM (n,e,u) = ('+str(lookN[i])+','+str(lookE[i])+','+str(lookU[i])+')'
print 'VLSM LOS = '+str(los_out[i])

#savetxt(out,c_[lon_out,lat_out,los_out,lookN,lookE,lookU],fmt='%14.6f\t%14.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f')