from matplotlib import pyplot as plt
from numpy import genfromtxt,argmin,array,zeros,ones,where,linspace,r_
from matplotlib.ticker import MultipleLocator


g=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat')
sta=genfromtxt('/Users/dmelgar/Amatrice2016/GPS/Cosismico_30Oct2016_GPS_GdL_V1.dat',usecols=0,dtype='S')
insar=genfromtxt(u'/Users/dmelgar/Amatrice2016/InSAR/M6.6/T44_Italy_1027_1102/T44_Italy_1027_1102.lltnde')

#Parse GPS
lon_gps=g[:,1]
lat_gps=g[:,2]
north=g[:,6]/1000
east=g[:,4]/1000
up=g[:,8]/1000


#parse insar
lon_insar=insar[:,0]
lat_insar=insar[:,1]
los=insar[:,6]/1000
lookE=insar[:,3]
lookN=insar[:,4]
lookU=insar[:,5]

#Projection variables
projected_gps=9999*ones(len(lon_gps))
los_insar=9999*ones(len(lon_gps))

thresh=0.005
for k in range(len(lon_gps)):
    
    #Get distance from GPS to LOS points
    d=((lon_gps[k]-lon_insar)**2+(lat_gps[k]-lat_insar)**2)**0.5
    i=argmin(d)
    
    if d[i]<thresh:
        
        #Get los vector
        #unit_vector=array([lookE[i],lookN[i],lookU[i]])
        unit_vector=array([lookE[i],lookN[i]])
        #project
        #projected_gps[k]=unit_vector.dot(array([east[k],north[k],up[k]]))
        projected_gps[k]=unit_vector.dot(array([east[k],north[k]]))
        los_insar[k]=los[i]
        
plt.figure(figsize=(6,10))
ax=plt.subplot(211)
plt.quiver(r_[11.65,lon_gps],r_[43.72,lat_gps],r_[0.1,east],r_[0,north],scale=1.2)
plt.annotate(s='10cm',xy=(11.65,43.82))
i=where(projected_gps<-0.2)[0]
#for k, txt in enumerate(sta[i]):
#    ax.annotate(txt, (lon_gps[i[k]],lat_gps[i[k]]))
ax.annotate('RIFP', (lon_gps[74],lat_gps[74]))
print (lon_gps[74],lat_gps[74])
ax.annotate('ARQT', (lon_gps[9],lat_gps[9]))
print (lon_gps[9],lat_gps[9])
#i=where(up<0)[0]
#j=where(up>=0)[0]
#plt.quiver(lon_gps[j],lat_gps[j],zeros(len(up[j])),up[j],scale=0.01,color='b')
#plt.quiver(lon_gps[i],lat_gps[i],zeros(len(up[i])),up[i],scale=0.01,color='r')

ax=plt.subplot(212)
i=where(projected_gps<9999)[0]
x=linspace(-0.4,0.15)
y=x
plt.plot(x,y,lw=2,c='k')
plt.scatter(projected_gps[i],los_insar[i],marker='s',s=30,lw=0.2,c='#0080FF')
i=where(projected_gps<-0.2)[0]
#for k, txt in enumerate(sta[i]):
#    ax.annotate(txt, (projected_gps[i[k]]+0.02,los_insar[i[k]]))
ax.annotate('RIFP', (projected_gps[74]+0.02,los_insar[74]))
ax.annotate('ARQT', (projected_gps[9]+0.02,los_insar[9]))
plt.xlim([-0.4,0.15])
plt.ylim([-0.4,0.15])
xmajorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
plt.ylabel('InSAR LOS (m)')
plt.xlabel('Projected GPS (m)')
plt.subplots_adjust(left=0.2,right=0.97,top=0.99,bottom=0.1)
plt.show()
