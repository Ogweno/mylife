from matplotlib import pyplot as plt
from numpy import genfromtxt,argmin,array,zeros,ones,where,linspace,r_
from matplotlib.ticker import MultipleLocator


g26=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat')
sta26=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_26Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
g30=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat')
sta30=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
insar_T44=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M6.6/T44_Italy_1027_1102/T44_Italy_1027_1102.lltnde')
insar_T22=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T22_Italy/T22_Italy.lltnde')
insar_T117=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M5.8_M6.6/T117_Italy/T117_Italy.lltnde')

#Parse GPS
gps=['VETT','CAMP','MSAN','ARQT','ASCC','ASC1','GNAL','MTER','MTTO','ACCU','GUMA','FIAB','GINE','CESI']

lon_gps=zeros(len(gps))
lat_gps=zeros(len(gps))

north26=zeros(len(gps))
east26=zeros(len(gps))
up26=zeros(len(gps))

north30=zeros(len(gps))
east30=zeros(len(gps))
up30=zeros(len(gps))

north2630=zeros(len(gps))
east2630=zeros(len(gps))
up2630=zeros(len(gps))

for k in range(len(gps)):
    
    i=where(sta26==gps[k])
    lon_gps[k]=g26[i,1]
    lat_gps[k]=g26[i,2]
    north26[k]=g26[i,6]/1000
    east26[k]=g26[i,4]/1000
    up26[k]=g26[i,8]/1000
    
    i=where(sta30==gps[k])
    north30[k]=g30[i,6]/1000
    east30[k]=g30[i,4]/1000
    up30[k]=g30[i,8]/1000
    
north2630=north26+north30
east2630=east26+east30
up2630=up26+up30

#parse insar
lon_insar_T44=insar_T44[:,0]
lat_insar_T44=insar_T44[:,1]
los_T44=insar_T44[:,6]/1000
lookE_T44=insar_T44[:,3]
lookN_T44=insar_T44[:,4]
lookU_T44=insar_T44[:,5]

lon_insar_T22=insar_T22[:,0]
lat_insar_T22=insar_T22[:,1]
los_T22=insar_T22[:,6]/1000
lookE_T22=insar_T22[:,3]
lookN_T22=insar_T22[:,4]
lookU_T22=insar_T22[:,5]

lon_insar_T117=insar_T117[:,0]
lat_insar_T117=insar_T117[:,1]
los_T117=insar_T117[:,6]/1000
lookE_T117=insar_T117[:,3]
lookN_T117=insar_T117[:,4]
lookU_T117=insar_T117[:,5]

#Projection variables
projected_gps_T44=9999*ones(len(lon_gps))
los_insar_T44=9999*ones(len(lon_gps))
projected_gps_T22=9999*ones(len(lon_gps))
los_insar_T22=9999*ones(len(lon_gps))
projected_gps_T117=9999*ones(len(lon_gps))
los_insar_T117=9999*ones(len(lon_gps))

thresh=0.005
for k in range(len(lon_gps)):
    
    #Get distance from GPS to LOS points
    d=((lon_gps[k]-lon_insar_T44)**2+(lat_gps[k]-lat_insar_T44)**2)**0.5
    i=argmin(d)
    
    if d[i]<thresh:
        
        #Get los vector
        unit_vector=array([lookE_T44[i],lookN_T44[i],lookU_T44[i]])
        #project
        projected_gps_T44[k]=unit_vector.dot(array([east30[k],north30[k],up30[k]]))
        los_insar_T44[k]=los_T44[i]
        
        
for k in range(len(lon_gps)):
    
    #Get distance from GPS to LOS points
    d=((lon_gps[k]-lon_insar_T22)**2+(lat_gps[k]-lat_insar_T22)**2)**0.5
    i=argmin(d)
    
    if d[i]<thresh:
        
        #Get los vector
        unit_vector=array([lookE_T22[i],lookN_T22[i],lookU_T22[i]])
        #project
        projected_gps_T22[k]=unit_vector.dot(array([east2630[k],north2630[k],up2630[k]]))
        los_insar_T22[k]=los_T22[i]
        
        
for k in range(len(lon_gps)):
    
    #Get distance from GPS to LOS points
    d=((lon_gps[k]-lon_insar_T117)**2+(lat_gps[k]-lat_insar_T117)**2)**0.5
    i=argmin(d)
    
    if d[i]<thresh:
        
        #Get los vector
        unit_vector=array([lookE_T117[i],lookN_T117[i],lookU_T117[i]])
        #project
        projected_gps_T117[k]=unit_vector.dot(array([east2630[k],north2630[k],up2630[k]]))
        los_insar_T117[k]=los_T117[i]
        
        
plt.figure(figsize=(10,6))
ax=plt.subplot(121)
plt.quiver(r_[13,lon_gps],r_[43.12,lat_gps],r_[0.1,east26],r_[0,north26],scale=0.8,color='red')
plt.quiver(r_[13,lon_gps],r_[43.12,lat_gps],r_[0.1,east30],r_[0,north30],scale=0.8,color='black')
plt.legend(['Oct 26','Oct 30'])
plt.annotate(s='10cm',xy=(12.96,43.15))
for k, txt in enumerate(gps):
    ax.annotate(txt, (lon_gps[k],lat_gps[k]))
    
ax=plt.subplot(122)
z=zeros(len(north26))
plt.quiver(r_[13,lon_gps],r_[43.10,lat_gps],r_[0,z],r_[0.1,up26],scale=0.8,color='red')
plt.quiver(r_[13,lon_gps],r_[43.10,lat_gps],r_[0,z],r_[0.1,up30],scale=0.8,color='black')
plt.legend(['Oct 26','Oct 30'])
plt.annotate(s='10cm',xy=(13.03,43.15))
for k, txt in enumerate(gps):
    ax.annotate(txt, (lon_gps[k],lat_gps[k]))

### Plot projections
plt.figure(figsize=(14,4.8))

ax=plt.subplot(131)
x=linspace(-0.4,0.4)
y=x
yplus=y+0.05
yminus=y-0.05
plt.plot(x,y,lw=2,c='k')
plt.plot(x,yplus,lw=1,c='k')
plt.plot(x,yminus,lw=1,c='k')
plt.scatter(projected_gps_T44,los_insar_T44,marker='s',s=30,lw=0.2,c='#0080FF')
for k, txt in enumerate(gps):
    ax.annotate(txt, (projected_gps_T44[k]+0.02,los_insar_T44[k]))
plt.xlim([-0.4,0.3])
plt.ylim([-0.3,0.2])
xmajorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
plt.ylabel('InSAR LOS T44 (m)')
plt.xlabel('Projected GPS Oct. 30th(m)')

ax=plt.subplot(132)
plt.plot(x,y,lw=2,c='k')
plt.plot(x,yplus,lw=1,c='k')
plt.plot(x,yminus,lw=1,c='k')
plt.scatter(projected_gps_T22,los_insar_T22,marker='s',s=30,lw=0.2,c='#FF8C00')
for k, txt in enumerate(gps):
    ax.annotate(txt, (projected_gps_T22[k]+0.02,los_insar_T22[k]))
plt.xlim([-0.4,0.3])
plt.ylim([-0.3,0.2])
xmajorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
plt.ylabel('InSAR LOS T22 (m)')
plt.xlabel('Projected GPS Oct. 26th+30th(m)')

ax=plt.subplot(133)
plt.plot(x,y,lw=2,c='k')
plt.plot(x,yplus,lw=1,c='k')
plt.plot(x,yminus,lw=1,c='k')
plt.scatter(projected_gps_T117,los_insar_T117,marker='s',s=30,lw=0.2,c='#228B22')
for k, txt in enumerate(gps):
    ax.annotate(txt, (projected_gps_T117[k]+0.02,los_insar_T117[k]))
plt.xlim([-0.4,0.3])
plt.ylim([-0.3,0.2])
xmajorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
plt.ylabel('InSAR LOS T117 (m)')
plt.xlabel('Projected GPS Oct. 26th+30th(m)')


plt.subplots_adjust(left=0.07,right=0.97,top=0.99,bottom=0.1,hspace=0.2)
plt.show()
