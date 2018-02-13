from numpy import genfromtxt,zeros,ones,savetxt,c_,array,sqrt,arange
from pyproj import Geod
from string import rjust

f=genfromtxt('/Users/dmelgar/Slip_inv/Slip_pulse/data/model_info/M8.fault')

vr=2.8
hypo=array([-90.,15.,16.])
 
rise_time=4.0
 
ss=zeros(len(f))
ds=ones(len(f))*2.8

#distance from subfault to hypo
p=Geod(ellps='WGS84')
az,baz,d=p.inv(f[:,1],f[:,2],ones(len(f))*hypo[0],ones(len(f))*hypo[1])
d=d/1000.
#add teh depth difference
d=d+sqrt((f[:,3]-hypo[2])**2)
t=d/vr

out=c_[f[:,0:7],rise_time*ones(len(f)),ss,ds,f[:,8:10],t,45e9*ones(len(f))]
fout='/Users/dmelgar/Slip_inv/Slip_pulse/output/forward_models/M8_4s.rupt'
savetxt(fout,out,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.1f\t%.6f\t%.6f\t%.1f\t%.1f\t%.1f\t%.1e')




##now make station profiles
#
##strike perpendicular
#lon1=arange(hypo[0]-0.5,hypo[0]+2,0.05)
#lat1=ones(len(lon1))*hypo[1]
#fout='/Users/dmelgar/Slip_inv/Slip_pulse/data/station_info/profiles.sta'
#
#f=open(fout,'w')
#for k in range(len(lon1)):
#    line='%s\t%.4f\t%.4f\n' % ('PE'+rjust(str(k),2,'0'),lon1[k],lat1[k])
#    f.write(line)
#    
##strike parallel
#lat2=arange(hypo[1],hypo[1]+2,0.05)
#lon2=ones(len(lat2))*hypo[0]
#for k in range(len(lon2)):
#    line='%s\t%.4f\t%.4f\n' % ('PA'+rjust(str(k),2,'0'),lon2[k],lat2[k])
#    f.write(line)
#f.close()