from numpy import genfromtxt,cos,deg2rad,ones,c_,arange,zeros,savetxt,r_
from pyproj import Geod
from mudpy import forward

P=Geod(ellps='WGS84')

col1=genfromtxt('/Users/dmelgar/Amatrice2016/Colfiorito1997/col1.fsp')
dx=2.5*1000
dy=2.5*1000
strike=152
dip=46

dist=dy*cos(deg2rad(dip))/2.
deltaz=dy*cos(deg2rad(dip))/2000.

#Get prpojected locations
lon1,lat1,az=P.fwd(col1[:,1],col1[:,0],ones(len(col1))*(strike+90),ones(len(col1))*dist)
z1=col1[:,4]+deltaz

#Put in rupture type array
# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
I=ones(len(col1))
Z=zeros(len(col1))
out1=c_[arange(1,len(z1)+1),lon1,lat1,z1,I*strike,I*dip,I,I,Z,-col1[:,5],I*dx,I*dy,Z,Z]



col2=genfromtxt('/Users/dmelgar/Amatrice2016/Colfiorito1997/col2.fsp')
dx=2.5*1000
dy=2.5*1000
strike=144
dip=42

dist=dy*cos(deg2rad(dip))/2.
deltaz=dy*cos(deg2rad(dip))/2000.

#Get prpojected locations
lon2,lat2,az=P.fwd(col2[:,1],col2[:,0],ones(len(col2))*(strike+90),ones(len(col2))*dist)
z2=col2[:,4]+deltaz

#Put in rupture type array
# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
I=ones(len(col2))
Z=zeros(len(col2))
out2=c_[arange(len(z1)+1,len(z1)+len(z2)+1),lon2,lat2,z2,I*strike,I*dip,I,I,Z,-col2[:,5],I*dx,I*dy,Z,Z]




col3=genfromtxt('/Users/dmelgar/Amatrice2016/Colfiorito1997/col3.fsp')
dx=1.5*1000
dy=1.5*1000
strike=135
dip=45

dist=dy*cos(deg2rad(dip))/2.
deltaz=dy*cos(deg2rad(dip))/2000.

#Get prpojected locations
lon3,lat3,az=P.fwd(col3[:,1],col3[:,0],ones(len(col3))*(strike+90),ones(len(col3))*dist)
z3=col3[:,4]+deltaz

#Put in rupture type array
# No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
I=ones(len(col3))
Z=zeros(len(col3))
out3=c_[arange(len(z1)+len(z2)+1,len(z1)+len(z2)+len(z3)+1),lon3,lat3,z3,I*strike,I*dip,I,I,Z,-col3[:,5],I*dx,I*dy,Z,Z]



out=r_[out1,out2,out3]
fout=u'/Users/dmelgar/Amatrice2016/Colfiorito1997/colfiorito.rupt'
fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
savetxt(fout,out,fmt=fmtout)

forward.inv2coulomb(u'/Users/dmelgar/Amatrice2016/Colfiorito1997/colfiorito.rupt',[13.3,42.6],'/Users/dmelgar/Amatrice2016/Colfiorito1997/colfiorito.coul')