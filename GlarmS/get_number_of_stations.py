from numpy import genfromtxt,linspace,ones,sin,cos,deg2rad,zeros,argmin,zeros,where,c_,savetxt,arange,meshgrid,savetxt
from gmsh_tools import llz2xyz
from pyproj import Geod
from scipy.interpolate import interp2d



#output where?
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_glarms.txt' #Panel 1
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_pbomerge.txt' #Panel 2
fout=u'/Users/dmelgar/GlarmS/Sensitivity/stations_over_2.0cm_M6.5_continuos.txt' #Panel 3

#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_glarms.txt' #Panel 1
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_pbomerge.txt' #Panel 2
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_continuos.txt' #Panel 3

#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.0_glarms.txt' #Panel 1
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.0_pbomerge.txt' #Panel 2
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.0_continuos.txt' #Panel 3

#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.5_glarms.txt' #Panel 1
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.5_pbomerge.txt' #Panel 2
#fout=u'/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M7.5_continuos.txt' #Panel 3

#regualr cascadia
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_20km.txt' 
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_30km.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_40km.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_50km.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_75km.txt'
#fout='/Users/dmelgar/GlarmS/Sensitivity/cascadia_stations_over_2.0cm_M6.5_100km.txt'



degkm=111.0

#read stations
#stations=genfromtxt('/Users/dmelgar/GlarmS/GlarmS_west_coast_stations.txt',usecols=[1,2]) #panel 1
#stations=genfromtxt('/Users/dmelgar/GlarmS/Glarm_pbo_merge.txt',usecols=[1,2])  #Panel 2
#stations=genfromtxt('/Users/dmelgar/GlarmS/stations_long_cascadia.txt') #Panel 3
stations=genfromtxt('/Users/dmelgar/GlarmS/stations_long.txt') #Panel 3

#read deformation
deformation=genfromtxt('/Users/dmelgar/GlarmS/Sensitivity/CA_M6.5_2.0cm.txt')
#deformation=genfromtxt('/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M6.5_2.0cm.txt')
#deformation=genfromtxt('/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M7.0_2.0cm.txt')
#deformation=genfromtxt('/Users/dmelgar/GlarmS/Sensitivity/Cascadia_M7.5_2.0cm.txt')


#Make grid of points I want to test
#Cascadia params
#nptsx=45
#nptsy=50
#lon_test=linspace(-128.5,-118,nptsx)
#lat_test=linspace(39,51,nptsy)
#uh_thres=0.02
#dist_thresh=2.0

#CA params
nptsx=40
nptsy=40
lon_test=linspace(-125,-114,nptsx)
lat_test=linspace(32,43,nptsy)
uh_thres=0.02
dist_thresh=2.0


#Make regularly sampled stations grid
#delta_km=100.0
#delta_degs=delta_km/degkm
#lon_sta=arange(-128.5,-118,delta_degs)
#lat_sta=arange(39,51,delta_degs)
#[X,Y]=meshgrid(lon_sta,lat_sta)
#lon_sta=X.ravel()
#lat_sta=Y.ravel()
#stations=zeros((len(lon_sta),2))
#stations[:,0]=lon_sta
#stations[:,1]=lat_sta
#savetxt(u'/Users/dmelgar/code/GMT/glarms_sensitivity/cascadia_stations_regular.txt',stations,fmt='%.6f')



#get cornrs of grid
xmin=deformation[:,0].min()
xmax=deformation[:,0].max()
ymin=deformation[:,1].min()
ymax=deformation[:,1].max()

p=Geod(ellps='WGS84')
def get_xy(stations,origin_lon,origin_lat):
    '''
    This function converts station coordiantes to x,y assuming some origin
    '''
    lo=ones(len(stations))*origin_lon
    la=ones(len(stations))*origin_lat
    az,baz,dist=p.inv(lo,la,stations[:,0],stations[:,1])
    dist=dist/1000.
    xstation=dist*sin(deg2rad(az))
    ystation=dist*cos(deg2rad(az))
    return xstation,ystation


uhsta=zeros(len(stations))
k=0
lon_out=zeros(len(lon_test)*len(lat_test))
lat_out=zeros(len(lon_test)*len(lat_test))
stations_over_threshold=zeros(len(lon_test)*len(lat_test))
print "I'm counting bro, hang on..."

for klon in range(len(lon_test)):
    for klat in range(len(lat_test)):
        
        if k%100 == 0:
            print '  %d of %d' % (k,len(lon_test)*len(lat_test))
        #Assume urrent test point is center of cartesian deformation grid get (x,y) of all sites
        origin_lon=lon_test[klon]
        origin_lat=lat_test[klat]
        
        #Convert to cartesian
        xsta,ysta=get_xy(stations,origin_lon,origin_lat)
        
        #keep only stations inside the grid
        i=where((xsta>xmin) & (xsta<xmax) & (ysta>ymin) & (ysta<ymax))[0]
        if len(i)>0:
            xsta=xsta[i]
            ysta=ysta[i]
            
            #interpolate horizontal displacement filed to station coordinates
            uh=(deformation[:,2]**2+deformation[:,3]**2)**0.5
            
            #For eeach station find closest grid point and make that be its uh
            for ksta in range(len(xsta)):
                dist=((deformation[:,0]-xsta[ksta])**2+(deformation[:,1]-ysta[ksta])**2)**0.5
                imin=argmin(dist)
                if dist[imin]>dist_thresh:
                    uhsta[i[ksta]]=0
                else:
                    uhsta[i[ksta]]=uh[imin]
            
            #count how mnay stations made the cut
            iuh=where(uhsta>uh_thres)[0]
            stations_over_threshold[k]=len(iuh)
        else:
            stations_over_threshold[k]=0
        lon_out[k]=origin_lon
        lat_out[k]=origin_lat
        k+=1

out=c_[lon_out,lat_out,stations_over_threshold]
savetxt(fout,out,fmt='%.6f\t%.6f\t%d')