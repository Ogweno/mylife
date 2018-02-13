from numpy import genfromtxt,where,linspace,savetxt,c_
from scipy.interpolate import interp1d
from string import rjust

#file made with:
# gmt pscoast -R-98/-90/11/17 -J -Di -A1000 -M -W > coast.txt

offset=6*15./3600

coast=genfromtxt('/Users/dmelgar/code/GMT/tehuantepec/coast.txt')


i=where(coast[:,1]<17)[0]
coast=coast[i,:]

i=where((coast[:,0]>-97)&(coast[:,0]<-91))[0]
coast=coast[i,:]

coast[:,1]=coast[:,1]-offset

#Now interpolate to every 0.05 degrees or wahtever
xi=linspace(coast[0,0],coast[-1,0],500)
f=interp1d(coast[:,0],coast[:,1])
yi=f(xi)

#Save for grdtrack
savetxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points.xy',c_[xi,yi],fmt='%.6f')

#Get grdtrack results
xyz=genfromtxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points.xyz')
i=where(xyz[:,2]<-5)[0]
xyz=xyz[i,:]
savetxt(u'/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points_filtered.xyz',xyz,fmt='%.4f')

#print gauge thingy to screen
f=open('/Users/dmelgar/DEMs/SRTM15/tehuant_coast_points.gauge','w')
for k in range(len(xyz)):
    num=rjust(str(k),3,'0')
    f.write('rundata.gaugedata.gauges.append([2'+num+','+str(xyz[k,0])+','+str(xyz[k,1])+',0., 1.e10])\n')
f.close()