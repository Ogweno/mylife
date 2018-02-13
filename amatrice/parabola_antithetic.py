from numpy import array,linspace,genfromtxt,r_,ones,where,c_,zeros
from scipy.linalg import inv
from matplotlib import pyplot as plt
from pyproj import Geod
import gmsh_tools

#load surface trace
surface_trace=genfromtxt(u'/Users/dmelgar/Amatrice2016/3D_fault/_listric/surface_trace_antithetic.txt')
gmsh_out=u'/Users/dmelgar/Amatrice2016/3D_fault/_listric/antithetic_listric.gmsh'
max_depth=-8.3
azimuth=155-90

#Load cross etion info
CC=genfromtxt('/Users/dmelgar/Amatrice2016/plots/CC.txt')
CCvettore=genfromtxt('/Users/dmelgar/Amatrice2016/plots/CCvettore.txt')
CCdetach=genfromtxt('/Users/dmelgar/Amatrice2016/plots/CCdetach.txt')
EE=genfromtxt('/Users/dmelgar/Amatrice2016/plots/EE.txt')
EEvettore=genfromtxt('/Users/dmelgar/Amatrice2016/plots/EEvettore.txt')
EEdetach=genfromtxt('/Users/dmelgar/Amatrice2016/plots/EEdetach.txt')
JJ=genfromtxt('/Users/dmelgar/Amatrice2016/plots/JJ.txt')
JJvettore=genfromtxt('/Users/dmelgar/Amatrice2016/plots/JJvettore.txt')
JJdetach=genfromtxt('/Users/dmelgar/Amatrice2016/plots/JJdetach.txt')

#Profiles
profiles=array([[ 13.01567151,  42.54648126,  13.51355551,  42.71660122],
       [ 12.98991192,  42.58726628,  13.48812229,  42.75738348],
       [ 12.96411872,  42.62804522,  13.4626561 ,  42.79815965],
       [ 12.93829183,  42.66881804,  13.43715684,  42.83892972],
       [ 12.91243116,  42.70958475,  13.41162445,  42.87969366],
       [ 12.88653663,  42.75034532,  13.38605883,  42.92045146],
       [ 12.86060814,  42.79109973,  13.36045989,  42.96120311],
       [ 12.83464561,  42.83184798,  13.33482755,  43.00194859],
       [ 12.80864894,  42.87259004,  13.30916172,  43.04268788],
       [ 12.78261806,  42.9133259 ,  13.28346231,  43.08342096]])
AAp=array([profiles[0,[0,2]],profiles[0,[1,3]]]).T
BBp=array([profiles[1,[0,2]],profiles[1,[1,3]]]).T
CCp=array([profiles[2,[0,2]],profiles[2,[1,3]]]).T
DDp=array([profiles[3,[0,2]],profiles[3,[1,3]]]).T
EEp=array([profiles[4,[0,2]],profiles[4,[1,3]]]).T
FFp=array([profiles[5,[0,2]],profiles[5,[1,3]]]).T
GGp=array([profiles[6,[0,2]],profiles[6,[1,3]]]).T
HHp=array([profiles[7,[0,2]],profiles[7,[1,3]]]).T
IIp=array([profiles[8,[0,2]],profiles[8,[1,3]]]).T
JJp=array([profiles[9,[0,2]],profiles[9,[1,3]]]).T


#EE is the main player define the 3 points x is along profile dist, y is depth
p1=[16.5,0]
p2=[18.6,-3.6]
p3=[21,-6]
p4=[23.5,-7]
p5=[27,-10]
x=[p1[0],p2[0],p3[0],p4[0],p5[0]]
y=[p1[1],p2[1],p3[1],p4[1],p5[1]]
xi=linspace(15,30,20)
G=array([[p1[0]**2,p1[0],1],[p2[0]**2,p2[0],1],[p3[0]**2,p3[0],1],[p4[0]**2,p4[0],1],[p5[0]**2,p5[0],1]])
d=array([[p1[1]],[p2[1]],[p3[1]],[p4[1]],[p5[1]]])
m=inv(G.T.dot(G)).dot(G.T).dot(d)
yi=m[0]*xi**2+m[1]*xi+m[2]
#Remove positive vlues of depth
i=where(yi<0)[0]
xi=xi[i]
yi=yi[i]
#Add zero
yi=r_[0,yi]
x0=(-m[1]-(m[1]**2-4*(m[0]*m[2]))**0.5)/(2*m[0])
xi=r_[x0,xi]

#MOVE TO CC PROFILE
delta=xi[-1]-31
xiCC=xi-delta

#MOVE TO JJ PROFILE
delta=xi[-1]-31.5
xiJJ=xi-delta

plt.figure()
plt.scatter(EE[:,0],EE[:,1],s=5,c='k')
plt.scatter(x,y,s=50)
plt.plot(xi,yi,lw=2)
plt.scatter(EEvettore[:,0],EEvettore[:,1],marker='s',s=60,c='r')
plt.scatter(EEdetach[:,0],EEdetach[:,1],marker='s',s=60,c='y')
plt.ylim([-15,0])

plt.figure()
plt.scatter(JJ[:,0],JJ[:,1],s=8,lw=0,c='magenta')
plt.plot(xiJJ,yi,lw=2)
plt.scatter(JJvettore[:,0],JJvettore[:,1],marker='s',s=60,c='r')
plt.scatter(JJdetach[:,0],JJdetach[:,1],marker='s',s=60,c='y')
plt.ylim([-15,0])

plt.figure()
plt.scatter(CC[:,0],CC[:,1],s=8,lw=0,c='blue')
plt.plot(xiCC,yi,lw=2,c='k')
plt.scatter(CCvettore[:,0],CCvettore[:,1],marker='s',s=60,c='r')
plt.scatter(CCdetach[:,0],CCdetach[:,1],marker='s',s=60,c='y')
plt.ylim([-15,0])
plt.show()


# Ok now get the parabolas at the edges of the surface trace
p=Geod(ellps='WGS84')
x_parabola=xi-xi[0]
y_parabola=yi
p0=surface_trace[0,:]
lon_parabola0,lat_parabola0,baz=p.fwd(p0[0]*ones(len(x_parabola)),p0[1]*ones(len(x_parabola)),azimuth*ones(len(x_parabola)),x_parabola*1000)
x_parabola=xi-xi[0]
y_parabola=yi
p1=surface_trace[-1,:]
lon_parabola1,lat_parabola1,baz=p.fwd(p1[0]*ones(len(x_parabola)),p1[1]*ones(len(x_parabola)),azimuth*ones(len(x_parabola)),x_parabola*1000)

#Filter by depth
i=where(y_parabola>max_depth)[0]
lon_parabola0=lon_parabola0[i]
lon_parabola1=lon_parabola1[i]
lat_parabola0=lat_parabola0[i]
lat_parabola1=lat_parabola1[i]
y_parabola=y_parabola[i]



#Now let's go to GMSH format brah
#Make a single variable with surface trace and two parabolas
fault=r_[c_[surface_trace,zeros(len(surface_trace))],c_[lon_parabola0,lat_parabola0,y_parabola],c_[lon_parabola1,lat_parabola1,y_parabola]]



#Move surface trace to bottom
#delta_x=surface_trace[-1,0]-lon_parabola1[1]
#delta_y=surface_trace[-1,1]-lat_parabola1[1]
#bottom_trace=zeros((len(surface_trace),3))
#bottom_trace[:,0]=surface_trace[:,0]-delta_x
#bottom_trace[:,1]=surface_trace[:,1]-delta_y
#bottom_trace[:,2]=ones(len(surface_trace))*fault[:,2].min()
#
##Append to fault
#fault=r_[fault,bottom_trace]

#Make gmsh file
gmsh_tools.xyz2gmsh(gmsh_out,fault[:,0],fault[:,1],fault[:,2],coord_type='UTM',projection_zone='33T')