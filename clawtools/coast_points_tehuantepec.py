'''
Here goes, read in a coast line file, filter to desired region, check support SRTM grid,
see if it's on a wet cell, if not then search for nearest wet cell.

'''

from numpy import r_,where,meshgrid,savetxt,c_,diff,array,arange
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from string import rjust

#Coquimbo
offset_integer=5
offset_integer1=2
delta=15./3600
offset=offset_integer*delta #2 for max ampl.

#srtm_file=u'/Users/dmelgar/DEMs/SRTM15/tehuantepec15.grd'
srtm_file=u'/Users/dmelgar/DEMs/SRTM15/tehuantepec_coast_points.grd'




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
#ilon=range(0,len(x),6)
#y=y[ilon]
#z=z[ilon,:]
srtm_lon,srtm_lat=meshgrid(x,y)

#Go line by line in srtm grid, finds bathy and get offfset point

lat_out=[]
lon_out=[]

#Go north to south
for klat in range(srtm_lon.shape[1]):
    k=klat
    zslice=z[:,k]
    ibathy=where(zslice<0)[0]
    itopo=where(zslice>0)[0]
    zslice[ibathy]=-1
    zslice[itopo]=1
    left=where(zslice>0)[0]
    right=where(zslice<0)[0]
        
    lat_out.append(srtm_lat[right[-1]-offset_integer,k])
    lon_out.append(srtm_lon[0,k-offset_integer1])

#To array

lon_out=array(lon_out)
lat_out=array(lat_out)

i=where(lat_out<16.218)[0]
lon_out=lon_out[i]
lat_out=lat_out[i]


#decimate by 3
i=arange(0,len(lon_out),3)

lon_out=lon_out[i]
lat_out=lat_out[i]


for k in range(len(lon_out)):
    num=rjust(str(k),3,'0')
    print 'rundata.gaugedata.gauges.append([1'+num+','+str(lon_out[k])+','+str(lat_out[k])+',0., 1.e10])'

#savetxt(fout,c_[lon_out,lat_out],fmt='%.8f\t%.8f')
            
