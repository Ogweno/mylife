from scipy.ndimage.filters import gaussian_filter
from numpy import genfromtxt,where,unique,meshgrid,zeros,savetxt,c_,unique,arange,ones,sqrt,argmin
from matplotlib import pyplot as plt


tsunami_in=genfromtxt(u'/Users/dmelgar/Coquimbo2015/tsunami/dtopo/gps_sm_tg_insar_final.dtopo')
#tsunami_in=genfromtxt(u'/Users/dmelgar/Coquimbo2015/tsunami/dtopo/uz.dtopo')
i=where(tsunami_in[:,0]>0)[0]
tsunami=tsunami_in[i,1:]
variance=1
coast=genfromtxt('/Users/dmelgar/Coquimbo2015/GV/coast_points.txt')
#Make profiles
lon_profiles=arange(-73.5,-71.2,0.01)
lat_profile1=ones(len(lon_profiles))*-29.89
lat_profile2=ones(len(lon_profiles))*-30.47
lat_profile3=ones(len(lon_profiles))*-31.20
#Look along coast
ucoast=zeros(len(coast))
for k in range(len(coast)):
    dist=sqrt((coast[k,0]-tsunami[:,0])**2+(coast[k,1]-tsunami[:,1])**2)
    imin=argmin(dist)
    ucoast[k]=tsunami[imin,2]
savetxt('/Users/dmelgar/Coquimbo2015/GV/uz_yeshoriz_grid.txt',c_[coast,ucoast*100],fmt='%12.6f\t%12.6f\t%10.4f',header='#lon,lat,uz(cm)')
#Now for profiles
uprofile1=zeros(len(lon_profiles))
uprofile2=zeros(len(lon_profiles))
uprofile3=zeros(len(lon_profiles))
for k in range(len(lon_profiles)):
    dist=sqrt((lon_profiles[k]-tsunami[:,0])**2+(lat_profile1[k]-tsunami[:,1])**2)
    imin=argmin(dist)
    uprofile1[k]=tsunami[imin,2]
    
    dist=sqrt((lon_profiles[k]-tsunami[:,0])**2+(lat_profile2[k]-tsunami[:,1])**2)
    imin=argmin(dist)
    uprofile2[k]=tsunami[imin,2]
    
    dist=sqrt((lon_profiles[k]-tsunami[:,0])**2+(lat_profile3[k]-tsunami[:,1])**2)
    imin=argmin(dist)
    uprofile3[k]=tsunami[imin,2]
savetxt('/Users/dmelgar/Coquimbo2015/GV/uz_yeshoriz_perfiles.txt',c_[lon_profiles,lat_profile1,uprofile1*100,lon_profiles,lat_profile2,uprofile2*100,lon_profiles,lat_profile3,uprofile3*100],fmt='%12.6f\t%12.6f\t%10.4f\t%12.6f\t%12.6f\t%10.4f\t%12.6f\t%12.6f\t%10.4f',header='#lon,lat,uz(m),lon,lat,uz(m),lon,lat,uz(m)')

plt.figure()
plt.scatter(tsunami[:,0],tsunami[:,1],c=tsunami[:,2],lw=0,vmin=-2,vmax=2,cmap='seismic')
plt.plot(lon_profiles,lat_profile1,lon_profiles,lat_profile2,lon_profiles,lat_profile3)
plt.scatter(coast[:,0],coast[:,1],c=ucoast,s=150,vmin=-2,vmax=2,cmap='seismic')
plt.colorbar()

plt.figure()
plt.plot(lon_profiles,uprofile1,lon_profiles,uprofile2,lon_profiles,uprofile3)
plt.legend(['29.89S','30.47S','31.20S'])
plt.xlabel('Longitude')
plt.ylabel('Vertical deformation (m)')

plt.show()