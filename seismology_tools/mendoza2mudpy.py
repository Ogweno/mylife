from numpy import genfromtxt,deg2rad,sin,cos,zeros,ones,arange,c_,savetxt
import pyproj
from string import rjust

sta_file='/Users/dmelgar/Michoacan1985/tsunami/mesh_res_05.xy'
sta_file_out='/Users/dmelgar/Michoacan1985/tsunami/seafloor_05.sta'
in_file='/Users/dmelgar/Michoacan1985/slip/mainshock_mendoza.txt'
out_file='/Users/dmelgar/Michoacan1985/slip/mich_main.fault'
out_rupt_file='/Users/dmelgar/Michoacan1985/slip/mich_main.rupt'

Nstrike=12
Ndip=10
Dx=15
Dy=13.9
dip=14
strike_proj=300
dip_proj=30
strike=300
rake=90
lon_top_left=-101.7616
lat_top_left=17.2604
z_top_left=6.1135

s=genfromtxt(in_file)/100
s=s.ravel()

#Get slip
ds=s*sin(deg2rad(rake))
ss=s*cos(deg2rad(rake))

#Get new coords
g = pyproj.Geod(ellps='WGS84')
lo=zeros(len(s))
la=zeros(len(s))
z=zeros(len(s))
k=0
for kdip in range(Ndip):
    for kstrike in range(Nstrike):
        #Go along strike
        lo1,la1,ba=g.fwd(lon_top_left,lat_top_left,strike_proj,(kstrike*Dx+Dx/2)*1000)
        #Go down-dip strike
        lo[k],la[k],ba=g.fwd(lo1,la1,dip_proj,(kdip*Dy+Dy/2)*1000)
        delta_z=Dy*(kdip+0.5)*sin(deg2rad(14))
        print delta_z
        z[k]=z_top_left+delta_z
        k+=1

#Write .rupt
i=ones(len(ds))
out=c_[arange(len(ds))+1,lo,la,z,i*strike,i*dip,0.5*i,i,ss,ds,Dx*i*1000,Dy*i*1000,i]
savetxt(out_rupt_file,out,fmt='%i\t%10.6f\t%10.6f\t%10.6f\t%7.2f\t%7.2f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%9.4f\t%9.4f\t%7.4f')

# Write .fault
i=ones(len(ds))
out=c_[arange(len(s))+1,lo,la,z,i*strike,i*dip,0.5*i,i,Dx*i*1000,Dy*i*1000]
savetxt(out_file,out,fmt='%i\t%10.6f\t%10.6f\t%10.6f\t%7.2f\t%7.2f\t%7.4f\t%7.4f\t%9.4f\t%9.4f')


##Now the station file
#lonlat=genfromtxt('/Users/dmelgar/Michoacan1985/tsunami/mesh_res_05.xy')
#f=open(sta_file_out,'w')
#for k in range(len(lonlat)):
#    line='SF'+rjust(str(k),4,'0')+'\t%10.6f\t%10.6f\n' % (lonlat[k,0],lonlat[k,1])
#    f.write(line)
#f.close()