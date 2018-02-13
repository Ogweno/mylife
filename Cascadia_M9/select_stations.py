from numpy import genfromtxt,zeros,where,unique
from string import rjust

dist_in_degs=2.0
fout=u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl'
hypocenter=[-124.616004,45.863800,19.84]
lonlat=genfromtxt(u'/Users/dmelgar/Cascadia_M9/eew_uw_prod1_channels.txt',usecols=[4,5])
sta=genfromtxt(u'/Users/dmelgar/Cascadia_M9/eew_uw_prod1_channels.txt',usecols=[1],dtype='S')

#Get stations based on distance
##Calcualte distance to all sites
#dist=((hypocenter[0]-lonlat[:,1])**2+(hypocenter[1]-lonlat[:,0])**2)**0.5
#i=where(dist<dist_in_degs)[0]
#lonlat=lonlat[i,:]
#sta=sta[i]
##get unique stations
#foo,j=unique(sta,return_index=True)
#lonlat=lonlat[j,:]
#sta=sta[j]

#get stations based on a box
lon=lonlat[:,1]
lat=lonlat[:,0]
i=where((lon>-126) & (lon<-123) & (lat>42.5) & (lat<46.5))[0]
#get unique stations
sta=sta[i]
lonlat=lonlat[i,:]
foo,j=unique(sta,return_index=True)
lonlat=lonlat[j,:]
sta=sta[j]




#writ eoutput
f=open(fout,'w')
for k in range(len(lonlat)):
    line='%.6f\t%.6f\t%s\t720\n' % (lonlat[k,1],lonlat[k,0],rjust(sta[k],4,' '))
    f.write(line)
f.close()