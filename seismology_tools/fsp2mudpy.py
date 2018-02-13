from numpy import genfromtxt,deg2rad,sin,cos,zeros,ones,arange,c_,savetxt
import pyproj
from string import rjust

sta_file='/Users/dmelgar/Michoacan1985/tsunami/mesh_res_05.xy'
sta_file_out='/Users/dmelgar/Michoacan1985/tsunami/seafloor_05.sta'
in_file='/Users/dmelgar/Michoacan1985/slip/s1985ZIHUAT01MEND.fsp'
out_file='/Users/dmelgar/Michoacan1985/slip/mich_after.fault'
out_rupt_file='/Users/dmelgar/Michoacan1985/slip/mich_after.rupt'

Dx=7.5
Dy=7.5
dip=14
co_strike=30
strike=300
rake=100

s=genfromtxt(in_file)
lon=s[:,1]
lat=s[:,0]

#Get slip
ds=s[:,5]*sin(deg2rad(rake))
ss=s[:,5]*cos(deg2rad(rake))

#Get new coords
z=s[:,4]+(Dy/2)*sin(deg2rad(dip))

#Reckon with pyproj
g = pyproj.Geod(ellps='WGS84')
#Now reckon
lo=zeros(len(s))
la=zeros(len(s))
for k in range(len(s)):
        lo[k],la[k],ba=g.fwd(lon[k],lat[k],co_strike,Dx*1000/2)

#Write .rupt
i=ones(len(s))
out=c_[arange(len(s))+1,lo,la,z,i*strike,i*dip,0.5*i,i,ss,ds,Dx*i*1000,Dy*i*1000,i]
savetxt(out_rupt_file,out,fmt='%i\t%10.6f\t%10.6f\t%10.6f\t%7.2f\t%7.2f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%9.4f\t%9.4f\t%7.4f')

# Write .fault
i=ones(len(s))
out=c_[arange(len(s))+1,lo,la,z,i*strike,i*dip,0.5*i,i,Dx*i*1000,Dy*i*1000]
savetxt(out_file,out,fmt='%i\t%10.6f\t%10.6f\t%10.6f\t%7.2f\t%7.2f\t%7.4f\t%7.4f\t%9.4f\t%9.4f')


#Now the station file
lonlat=genfromtxt('/Users/dmelgar/Michoacan1985/tsunami/mesh_res_05.xy')
f=open(sta_file_out,'w')
for k in range(len(lonlat)):
    line='SF'+rjust(str(k),4,'0')+'\t%10.6f\t%10.6f\n' % (lonlat[k,0],lonlat[k,1])
    f.write(line)
f.close()