from numpy import savetxt,linspace,zeros,ones,arange,c_

fout='/Users/dmelgar/Amatrice2016/mt3d/big_volume.sources'
lon1=12.9
lon2=13.4
lat1=42.6
lat2=43.1
z1=1
z2=13

nx=30
ny=30
nz=13

x=linspace(lon1,lon2,ny)
y=linspace(lat1,lat2,nx)
z=linspace(z1,z2,nz)

lon=zeros(nx*ny*nz)
lat=zeros(nx*ny*nz)
depth=zeros(nx*ny*nz)

k=0
for kz in range(nz):
    for klat in range(ny):
        for klon in range(nx):
            lon[k]=x[klon]
            lat[k]=y[klat]
            depth[k]=z[kz]
            k+=1
            

savetxt(fout,c_[arange(1,len(lon)+1),lon,lat,depth],fmt='%d\t%.6f\t%.6f\t%.6f',header='No., lon, lat, z(km)')
        