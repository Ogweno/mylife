from pyproj import Geod
from numpy import cos,deg2rad,r_,array,c_,savetxt

P=Geod(ellps='WGS84')
#
##Colfiorito 1
#fout=u'/Users/dmelgar/code/GMT/Amatrice/colfiorito1_polygon.txt'
#lat=43.0225
#lon=12.8917
#strike=152
#dip=46
#length=7.5*1e3
#width=7.5*1e3
#
##
#P1=[lon,lat]
#lo,la,baz=P.fwd(lon,lat,strike,length)
#P2=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+90,width*cos(deg2rad(dip)))
#P3=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+180,length)
#P4=[lo,la]
#P5=[lon,lat]
#out=c_[array(P1),array(P2),array(P3),array(P4),array(P5)].T
#savetxt(fout,out,fmt='%.4f')
#
#
#
##Colfiorito 2
#fout=u'/Users/dmelgar/code/GMT/Amatrice/colfiorito2_polygon.txt'
#lat=43.0305
#lon=12.8622
#strike=144
#dip=42
#length=12.5*1e3
#width=7.5*1e3
#
##
#P1=[lon,lat]
#lo,la,baz=P.fwd(lon,lat,strike,length)
#P2=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+90,width*cos(deg2rad(dip)))
#P3=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+180,length)
#P4=[lo,la]
#P5=[lon,lat]
#out=c_[array(P1),array(P2),array(P3),array(P4),array(P5)].T
#savetxt(fout,out,fmt='%.4f')
#
#
#
#
##Colfiorito 3
#fout=u'/Users/dmelgar/code/GMT/Amatrice/colfiorito3_polygon.txt'
#lat=42.919
#lon=12.926
#strike=135
#dip=45
#length=9*1e3
#width=6*1e3
#
##
#P1=[lon,lat]
#lo,la,baz=P.fwd(lon,lat,strike,length)
#P2=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+90,width*cos(deg2rad(dip)))
#P3=[lo,la]
#lo,la,baz=P.fwd(lo,la,strike+180,length)
#P4=[lo,la]
#P5=[lon,lat]
#out=c_[array(P1),array(P2),array(P3),array(P4),array(P5)].T
#savetxt(fout,out,fmt='%.4f')
#



#Norica
fout=u'/Users/dmelgar/code/GMT/Amatrice/Norcia_east_polygon.txt'
lat=42.7255
lon=13.0403
strike=320
dip=70
length=9.8*1e3
width=12.5*1e3

#
P1=[lon,lat]
lo,la,baz=P.fwd(lon,lat,strike,length)
P2=[lo,la]
lo,la,baz=P.fwd(lo,la,strike+90,width*cos(deg2rad(dip)))
P3=[lo,la]
lo,la,baz=P.fwd(lo,la,strike+180,length)
P4=[lo,la]
P5=[lon,lat]
out=c_[array(P1),array(P2),array(P3),array(P4),array(P5)].T
savetxt(fout,out,fmt='%.4f')