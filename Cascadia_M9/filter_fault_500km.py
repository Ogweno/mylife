from numpy import genfromtxt,ones,where,savetxt
from pyproj import Geod

dlim=490.0 #in km

#read stations
sta_lonlat=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.stl',usecols=[0,1])

#Read fault model
fault=genfromtxt(u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000009.rupt')
fout=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/cascadia.planar.short.rupt'

#Calcualte distances from subfault to alls tations if distance?dlim give it a 0
p=Geod(ellps='WGS84')
keep=[]
for k in range(len(fault)):
    fault_lon=ones(len(sta_lonlat))*fault[k,1]
    fault_lat=ones(len(sta_lonlat))*fault[k,2]
    az,baz,dist=p.inv(sta_lonlat[:,0],sta_lonlat[:,1],fault_lon,fault_lat)
    i=where(dist/1000.>dlim)[0]
    if i.sum()<1: #one station is over threshold
        keep.append(k)
        
#Save
#1	-125.886749	 49.186325	 27.3840	 318.45	  10.61	 0.5	 9.32	 0.00	 4.99	   3451.09	   3451.09	117.58	4.663229e+10

savetxt(fout,fault[keep,:],fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e',header='#')
        