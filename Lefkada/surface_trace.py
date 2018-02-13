from numpy import genfromtxt,c_,ones,zeros,array,cross,unique,savetxt
from numpy.linalg import lstsq
 
f=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada/data/model_info/lefkada80.fault')
x=f[:,1]
y=f[:,2]
z=f[:,3]

#define 3 points
P=array([x[0],y[0],z[0]])
Q=array([x[1],y[1],z[1]])
R=array([x[-1],y[-1],z[-1]])
#2 vectors
PQ=Q-P
PR=R-P
#normal vector
n=cross(PQ,PR)
# First 3 parameters
A=n[0]
B=n[1]
C=n[2]
D=-P[0]*A-P[1]*B-P[2]*C
#Get equation of surface intersection
lon_surface=x[0:40]
m=-A/B
b=-D/B
lat_surface=m*lon_surface+b
#
savetxt('/Users/dmelgar/Slip_inv/Lefkada/data/model_info/lefkada80_surface_trace.txt',c_[lon_surface,lat_surface])