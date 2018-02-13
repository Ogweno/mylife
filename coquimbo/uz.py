from scipy.ndimage.filters import gaussian_filter
from numpy import genfromtxt,where,unique,meshgrid,zeros,savetxt,c_,unique,arange,ones
from scipy.interpolate import interp2d
from matplotlib import pyplot as plt


#tsunami=genfromtxt(u'/Users/dmelgar/Coquimbo2015/tsunami/dtopo/gps_sm_tg_insar_final.dtopo')
tsunami_in=genfromtxt(u'/Users/dmelgar/Coquimbo2015/tsunami/dtopo/uz.dtopo')
i=where(tsunami_in[:,0]>0)[0]
tsunami=tsunami_in[i,1:]
variance=1
coast=genfromtxt('/Users/dmelgar/Coquimbo2015/GV/coast_points.txt')
#Make profiles
lon_profiles=arange(-73.5,-71.2,0.01)
lat_profile1=ones(len(lon_profiles))*-29.89
lat_profile2=ones(len(lon_profiles))*-30.47
lat_profile3=ones(len(lon_profiles))*-31.20

[x,y]=meshgrid(unique(tsunami[:,0]),unique(tsunami[:,1]))
k=0
uz=zeros(x.shape)
for i in range(len(unique(tsunami[:,1]))):
    for j in range(len(unique(tsunami[:,0]))):
        uz[i,j]=tsunami[k,2]
        k+=1
f = interp2d(unique(tsunami[:,0]), unique(tsunami[:,1]), uz, kind='cubic')

#Test integrity
#j=arange(0,221328,10)
#u=zeros(len(j))
#for k in range(len(j)):
#    u[k]=f(tsunami[j[k],0],tsunami[j[k],1])

#uout=gaussian_filter(uz,variance)
ucoast=zeros(len(coast))
for k in range(len(coast)):
    ucoast[k]=f(coast[k,0],coast[k,1])
savetxt('/Users/dmelgar/Coquimbo2015/GV/uz_nohoriz.txt',c_[coast,ucoast*100],fmt='%12.6f\t%12.6f\t%10.4f',header='#lon,lat,uz(m)')
#Now for profiles
uprofile1=zeros(len(lon_profiles))
uprofile2=zeros(len(lon_profiles))
uprofile3=zeros(len(lon_profiles))
for k in range(len(lon_profiles)):
    uprofile1[k]=f(lon_profiles[k],lat_profile1[k])
    uprofile2[k]=f(lon_profiles[k],lat_profile2[k])
    uprofile3[k]=f(lon_profiles[k],lat_profile3[k])
savetxt('/Users/dmelgar/Coquimbo2015/GV/uz_nohoriz_perfiles.txt',c_[lon_profiles,lat_profile1,uprofile1*100,lon_profiles,lat_profile2,uprofile2*100,lon_profiles,lat_profile3,uprofile3*100],fmt='%12.6f\t%12.6f\t%10.4f\t%12.6f\t%12.6f\t%10.4f\t%12.6f\t%12.6f\t%10.4f',header='#lon,lat,uz(m),lon,lat,uz(m),lon,lat,uz(m)')

plt.figure()
plt.plot(lon_profiles,uprofile1,lon_profiles,uprofile2,lon_profiles,uprofile3)
plt.legend(['29.89S','30.47S','31.20S'])
plt.xlabel('Longitude')
plt.ylabel('Vertical deformation (m)')
plt.show()
