from numpy import genfromtxt,ones,r_,c_,tan,deg2rad,zeros,savetxt,flipud,expand_dims
import gmsh_tools
from pyproj import Geod

#gmsh_out='/Users/dmelgar/Amatrice2016/3D_fault/amatrice_long.gmsh'
#gmsh_out='/Users/dmelgar/Amatrice2016/3D_fault/amatrice_long.gmsh'
surface_trace=genfromtxt('/Users/dmelgar/Amatrice2016/3D_fault/Antithetic_fault_trace.txt')
gmsh_out='/Users/dmelgar/Amatrice2016/3D_fault/Antithetic_fault.gmsh'
horizontal_distance=12.0*ones(len(surface_trace)) #in km
dip=50

#Get average "strike" by getting azimuth between start and end poins
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(surface_trace[0,0],surface_trace[0,1],surface_trace[-1,0],surface_trace[-1,1])
proj_direction=(az+90)*ones(len(surface_trace))
lon_deep,lat_deep,baz=g.fwd(surface_trace[:,0],surface_trace[:,1],proj_direction,horizontal_distance*1000)
z_deep=-horizontal_distance*tan(deg2rad(dip))


#Make a single variable
fault=r_[c_[surface_trace,zeros(len(surface_trace))],c_[lon_deep,lat_deep,z_deep]]


#poly
npoints=len(fault)
poly=r_[fault[0:npoints/2,0:2],flipud(fault[npoints/2:-1,0:2]),expand_dims(fault[0,0:2],0)]
savetxt(u'/Users/dmelgar/Amatrice2016/3D_fault/polygon_antithetic.txt',poly,fmt='%.4f\t%.4f')



#Make gmsh file
gmsh_tools.xyz2gmsh(gmsh_out,fault[:,0],fault[:,1],fault[:,2],coord_type='UTM',projection_zone='33T')
