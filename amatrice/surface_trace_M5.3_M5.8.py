from numpy import genfromtxt,ones,deg2rad,tan
from pyproj import Geod
from matplotlib import pyplot as plt

deep_depth=15.0
surface_trace=genfromtxt('/Users/dmelgar/Amatrice2016/surface_traces/M5.3-M5.8.txt')
expected_depth=deep_depth*ones(len(surface_trace)) #in km
dip=50
horiz_deep=expected_depth/tan(deg2rad(dip))

#Get average "strike" by getting azimuth between start and end poins
g=Geod(ellps='WGS84')
strike,baz,dist=g.inv(surface_trace[0,0],surface_trace[0,1],surface_trace[-1,0],surface_trace[-1,1])
proj_direction=(strike+90)*ones(len(surface_trace))
lon_deep,lat_deep,baz=g.fwd(surface_trace[:,0],surface_trace[:,1],proj_direction,horiz_deep*1000)
print 'Strike is '+str(strike)

#Get middle of rectangle
lon_top_middle,lat_top_middle,baz=g.fwd(surface_trace[0,0],surface_trace[0,1],strike,dist/2)

#And go down dip
lon_middle,lat_middle,baz=g.fwd(lon_top_middle,lat_top_middle,strike+90,horiz_deep[0]*1000/2)

#Depth of the middle
print 'hypo is lon=%.4f, lat=%.4f, z=%.2f' % (lon_middle,lat_middle,deep_depth/2)


#
plt.scatter(surface_trace[:,0],surface_trace[:,1])
plt.scatter(lon_deep,lat_deep)
plt.scatter(lon_middle,lat_middle,c='r')
plt.axis('equal')
plt.show()