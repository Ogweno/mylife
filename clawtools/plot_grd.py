'''
Quick plot of any GMT grd file, this will only work for NETCDF4 files, 
i.e. if you use GMT5. If you are outdated and use NETCDF3 you can edit this
to use scipy.io.netcdf_file instead.

grdfile - path to file
'''
from netCDF4 import Dataset
from numpy import meshgrid,genfromtxt
import matplotlib.pyplot as plt
from matplotlib import cm

grdfile=u'/Users/dmelgar/DEMs/SRTM15/michoacan_srtm15.grd'
zlims=[-10,10]
cmap=cm.seismic
flip_lon=False

grd = Dataset(grdfile, 'r', format='NETCDF4')
try:
    x=grd.variables['x'][:]
    y=grd.variables['y'][:]
    z=grd.variables['z'][:]
except:
    x=grd.variables['lon'][:]
    y=grd.variables['lat'][:]
    z=grd.variables['z'][:]
if flip_lon:
    x=x-360
X,Y=meshgrid(x,y)
plt.figure()
plt.title(grdfile)
plt.pcolormesh(X,Y,z,vmin=zlims[0],vmax=zlims[1],cmap=cmap)
plt.colorbar()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()