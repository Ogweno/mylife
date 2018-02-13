from numpy import genfromtxt,ones,deg2rad,tan
from pyproj import Geod
from matplotlib import pyplot as plt

surface_trace=genfromtxt('/Users/dmelgar/Amatrice2016/3D_fault/surface_trace_fits_geol.txt')
expected_depth=8.0*ones(len(surface_trace)) #in km
dip=52
horiz_deep=expected_depth/tan(deg2rad(dip))

#Get average "strike" by getting azimuth between start and end poins
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(surface_trace[0,0],surface_trace[0,1],surface_trace[-1,0],surface_trace[-1,1])
proj_direction=(az+90)*ones(len(surface_trace))
lon_deep,lat_deep,baz=g.fwd(surface_trace[:,0],surface_trace[:,1],proj_direction,horiz_deep*1000)


#
plt.scatter(surface_trace[:,0],surface_trace[:,1])
plt.scatter(lon_deep,lat_deep)
plt.scatter(13.214403,42.719490,c='r')
plt.axis('equal')
plt.show()