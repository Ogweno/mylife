'''
Here goes, read in a coast line file, filter to desired region, check support SRTM grid,
see if it's on a wet cell, if not then search for nearest wet cell.

'''

from numpy import r_,where,meshgrid,savetxt,c_,diff,array,arange
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from string import rjust

#Coquimbo
offset_integer=2
delta=15./3600
offset=offset_integer*delta #2 for max ampl.

srtm_file=u'/Users/dmelgar/DEMs/SRTM15/kaikoura.grd'
fout=u'/Users/dmelgar/NewZealand2016/tsunami/kaikoura_coast_fg_text.txt'
fout2=u'/Users/dmelgar/NewZealand2016/tsunami/kaikoura_coast_fg_text.txt'


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
ilat=range(0,len(y),6)
y=y[ilat]
z=z[ilat,:]
srtm_lon,srtm_lat=meshgrid(x,y)

#Go line by line in srtm grid, finds bathy and get offfset point

lat_out=[]
lon_out=[]
#Just look for negative points
#for k in range(srtm_lat.shape[0]):
#    zslice=z[k,:]
#    ibathy=where(zslice<0)[0]
#    lon_out[k]=srtm_lon[k,ibathy[-1]-offset]
#Go left to right
for k in range(srtm_lat.shape[0]):
    zslice=z[k,:]
    ibathy=where(zslice<0)[0]
    itopo=where(zslice>0)[0]
    zslice[ibathy]=-1
    zslice[itopo]=1
    zdiff=diff(zslice)
    left=where(zdiff>0)[0]
    right=where(zdiff<0)[0]
    
    #for j in range(len(left)):
    #    #is point wet?
    #    if z[k,left[j]-offset_integer]<0:
    #        lat_out.append(srtm_lat[k][0])
    #        lon_out.append(srtm_lon[k,left[j]]-offset)
        
    for j in range(len(right)):
        #is point wet
        if z[k,right[j]+offset_integer]<0:
            lat_out.append(srtm_lat[k][0])
            lon_out.append(srtm_lon[k,right[j]]+offset)

#To array
lon_out=array(lon_out)
lat_out=array(lat_out)


#line filter
m=(-44.4857+39.5484)/(170.56-176.9)
b=-44.4857-m*170.56
i=where(lat_out<(lon_out*m+b))[0]

lon_out=lon_out[i]
lat_out=lat_out[i]


i=where((lon_out>171) & (lon_out<178) & (lat_out>-45) & (lat_out<-40))[0]
lon_out=lon_out[i]
lat_out=lat_out[i]

f=open(fout,'w')
for k in range(len(lon_out)):
    f.write('rundata.gaugedata.gauges.append([1'+rjust(str(k),3,'0')+','+str(lon_out[k])+','+str(lat_out[k])+',0., 1.e10])\n')
f.close()

savetxt(fout2,c_[lon_out,lat_out],fmt='%.8f\t%.8f')
            
