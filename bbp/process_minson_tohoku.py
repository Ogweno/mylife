from numpy import genfromtxt,where,cos,deg2rad,ones,sin,r_,expand_dims
from pyproj import Geod
from gmsh_tools import xyz2gmsh

gmsh_out='/Users/dmelgar/Tohoku2011/Minsons/rupt.gmsh'

dl=841**0.5

f=genfromtxt('/Users/dmelgar/Tohoku2011/Minsons/rupt.txt')
strike=194.
perp_strike=strike-90.
shallow_dip=3.
deep_dip=29.

p=Geod(ellps='WGS84')

#Get shallowest points
i=where(f[:,3]==8.07)[0]
shallow=f[i,:].copy()
dist=(dl/2.)*cos(deg2rad(shallow_dip))
shallow[:,2],shallow[:,1],foo=p.fwd(shallow[:,2],shallow[:,1],ones(len(shallow))*perp_strike,ones(len(shallow))*dist*1000)
shallow[:,3]=shallow[:,3]-(dl/2.)**sin(deg2rad(shallow_dip))
shallow=r_[shallow[:,:],expand_dims(shallow[-1,:],0),expand_dims(shallow[-1,:],0)]
shallow[-1,2],shallow[-1,1],foo=p.fwd(shallow[0,2],shallow[0,1],strike-180,dl*1000./2)
shallow[-2,2],shallow[-2,1],foo=p.fwd(shallow[-2,2],shallow[-2,1],strike,dl*1000./2)

#Get deepest points
i=where(f[:,3]==63.94)[0]
deep=f[i,:].copy()
dist=(dl/2.)*cos(deg2rad(deep_dip))
deep[:,2],deep[:,1],foo=p.fwd(deep[:,2],deep[:,1],ones(len(deep))*perp_strike+180,ones(len(deep))*dist*1000)
deep[:,3]=deep[:,3]+(dl/2.)**sin(deg2rad(deep_dip))
deep=r_[deep[:,:],expand_dims(deep[-1,:],0),expand_dims(deep[-1,:],0)]
deep[-1,2],deep[-1,1],foo=p.fwd(deep[0,2],deep[0,1],strike-180,dl*1000./2)
deep[-2,2],deep[-2,1],foo=p.fwd(deep[-2,2],deep[-2,1],strike,dl*1000./2)

#left
left=f[0:9,:].copy()
left[:,2],left[:,1],foo=p.fwd(left[:,2],left[:,1],ones(len(left))*strike-180,ones(len(left))*0.5*dl*1000)

#Right
right=f[207:216,:].copy()
right[:,2],right[:,1],foo=p.fwd(right[:,2],right[:,1],ones(len(right))*strike,ones(len(right))*0.5*dl*1000)

out=r_[f,shallow,deep,right,left]

#Gmshize

xyz2gmsh(gmsh_out,out[:,2],out[:,1],-out[:,3],coord_type='UTM',projection_zone='54S')