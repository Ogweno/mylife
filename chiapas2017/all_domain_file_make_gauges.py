from numpy import genfromtxt,where,meshgrid,savetxt,c_,arange
from scipy.interpolate import interp1d
from string import rjust

#file made with:
# gmt pscoast -R-98/-90/11/17 -J -Di -A1000 -M -W > coast.txt

dx=0.05
x=arange(-97,-91,0.1)
y=arange(13,17,0.1)

X,Y=meshgrid(x,y)
x=X.ravel()
y=Y.ravel()


#Save for grdtrack
savetxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points.xy',c_[x,y],fmt='%.6f')

#Get grdtrack results
xyz=genfromtxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points.xyz')
i=where(xyz[:,2]<-5)[0]
xyz=xyz[i,:]
savetxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points_inwater.xyz',xyz,fmt='%.4f')

#print gauge thingy to screen
f=open('/Users/dmelgar/DEMs/SRTM15/tehuant_domain_points.gauge','w')
for k in range(len(xyz)):
    num=rjust(str(k),4,'0')
    f.write('rundata.gaugedata.gauges.append([9'+num+','+str(xyz[k,0])+','+str(xyz[k,1])+',0., 1.e10])\n')
f.close()