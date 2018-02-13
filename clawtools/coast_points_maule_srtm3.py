'''
Here goes, read in a coast line file, filter to desired region, check support SRTM grid,
see if it's on a wet cell, if not then search for nearest wet cell.

'''

from numpy import r_,where,meshgrid,reshape,argsort,squeeze,linspace,savetxt,c_,diff,zeros,arange
from netCDF4 import Dataset
from scipy.interpolate import interp1d

##iquique_gps_long
#offset=5  #2 for max ampl.
#srtm_file=u'/Users/dmelgar/DEMs/SRTM30/iquique_claw.grd'
#fout=u'/Users/dmelgar/Iquique2014/tsunami/iquique_claw_30_times.txt'

#iquique_gps
#offset=3  #2 for max ampl.
#srtm_file=u'/Users/dmelgar/DEMs/iquique/iquique_3_final.grd'
#fout=u'/Users/dmelgar/Iquique2014/tsunami/iquique_claw_3_coast.txt'
#Maule GPS 30s
offset=3  #2 for max ampl.
srtm_file=u'/Users/dmelgar/DEMs/maule/maule_3_final.grd'
fout=u'/Users/dmelgar/Maule2010/tsunami/maule_3_coast_fg.txt'


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
ilat=range(0,len(y),10)
y=y[ilat]
z=z[ilat,:]
srtm_lon,srtm_lat=meshgrid(x,y)

#Go line by line in srtm grid, finds bathy and get offfset point

lat_out=y.copy()
lon_out=zeros(lat_out.shape)
#Just look for negative points
#for k in range(srtm_lat.shape[0]):
#    zslice=z[k,:]
#    ibathy=where(zslice<0)[0]
#    lon_out[k]=srtm_lon[k,ibathy[-1]-offset]
#Go left to right
for k in range(srtm_lat.shape[0]):
    zslice=z[k,:]
    ibathy=where(zslice<0)[0]
    ibathy_diff=diff(ibathy)
    ifirst=where(ibathy_diff>1)[0]#This is the las contiguos point from left to right, jesus
    if len(ifirst)==0: #no anomalies pick last negative
        ifirst=ibathy[-1]
    else: #Pick last contiguos
        ifirst=ifirst[0]
    lon_out[k]=srtm_lon[k,ifirst-offset]
#Decimate
idec=arange(0,len(lon_out),2)
lon_out=lon_out[idec]
lat_out=lat_out[idec]
    

savetxt(fout,c_[lon_out,lat_out],fmt='%.8f\t%.8f')
            
