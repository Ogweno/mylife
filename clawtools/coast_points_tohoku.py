'''
Here goes, read in a coast line file, filter to desired region, check support SRTM grid,
see if it's on a wet cell, if not then search for nearest wet cell.

'''

from numpy import r_,where,meshgrid,savetxt,c_,arange,array,argsort
from netCDF4 import Dataset


offset=2  #2 for max ampl.
#srtm_file=u'/Users/dmelgar/DEMs/tohoku/tohoku30.grd'
#fout=u'/Users/dmelgar/Tohoku2011/tsunami/tohoku_30_coast.fg'
dec=1
srtm_file=u'/Users/dmelgar/DEMs/tohoku/tohoku3_combined.grd'
fout=u'/Users/dmelgar/Tohoku2011/tsunami/tohoku_3_coast.fg'
dec=15
blend=True

#Read SRTM data
srtm = Dataset(srtm_file, 'r', format='NETCDF4')
try:
    x=srtm.variables['x'][:]
    y=srtm.variables['y'][:]
    z=srtm.variables['z'][:]
except:
    x=srtm.variables['lon'][:]
    y=srtm.variables['lat'][:]
    z=srtm.variables['z'][:]
ilat=range(0,len(y),dec)
y=y[ilat]
z=z[ilat,:]
srtm_lon,srtm_lat=meshgrid(x,y)

#Go line by line in srtm grid, finds bathy and get offfset point

lat_out=[]
lon_out=[]

for k in range(srtm_lat.shape[0]):
    print k
    zslice=z[k,:]
    #Brute force it, look from right to left
    i=where(zslice>0)[0]
    if len(i)>0:
        position=i[-1] #Last one from left to right
        if position+offset<srtm_lon.shape[1]:
            lon_out.append(srtm_lon[k,position+offset])
            lat_out.append(srtm_lat[k,0])
lat_out=array(lat_out)
lon_out=array(lon_out)
# Now go north south
lat_out2=[]
lon_out2=[]

for k in range(srtm_lon.shape[1]):
    print k
    zslice=z[:,k]
    #Brute force it, look from right to left
    i=where(zslice>0)[0]
    if len(i)>0:
        position=i[0] #First one from south to north
        if position-offset>0:
            lon_out2.append(srtm_lon[0,k])
            lat_out2.append(srtm_lat[position-offset,k])
lat_out2=array(lat_out2)
lon_out2=array(lon_out2)
#Blend
if blend==True:
    lon_out=r_[lon_out,lon_out2]
    lat_out=r_[lat_out,lat_out2]
    i=argsort(lat_out)
    lat_out=lat_out[i]
    lon_out=lon_out[i]
 
#Decimate=
lon_out=array(lon_out)
lat_out=array(lat_out)
idec=arange(0,len(lat_out),1)
lon_out=lon_out[idec]
lat_out=lat_out[idec]
#    

savetxt(fout,c_[lon_out,lat_out],fmt='%.8f\t%.8f')
            
