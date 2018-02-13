from numpy import genfromtxt,ones,where,savetxt,arange,r_,sort
from pyproj import Geod
from mudpy import forward
from matplotlib import pyplot as plt
from pyproj import Geod

#Source is 420x52
nrows=52
ncolumns=420
N=nrows*ncolumns

dlim=490.0 #in km

#Read fault model
fault=genfromtxt(u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.rupt')
# Need log for magnitude and hypocenter NOTHING ELSE
log=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.000006.log'
#Output rupt file
fout=u'/Users/dmelgar/FakeQuakes/Cascadia_M9/output/ruptures/planar.short.000006.rupt'


#crop columns
keep_columns=arange(100,328)
for k in range(len(keep_columns)):
    if k==0:
        i=arange(keep_columns[k],N,ncolumns)
    else:
        i=r_[i,arange(keep_columns[k],N,ncolumns)]
i=sort(i)
fault=fault[i,:]
fault[:,0]=arange(1,len(fault)+1)

#load stations
sta_lonlat=genfromtxt(u'/Users/dmelgar/Cascadia_M9/pnw.planar.0006.stl',usecols=[0,1])

#Find maximum subfault-station distance


#make plot
plt.figure()
plt.scatter(fault[:,1],fault[:,2],c=fault[:,9],lw=0,s=5)
plt.scatter(sta_lonlat[:,0],sta_lonlat[:,1],c='w')
plt.axis('equal')
plt.show()

#get longest distance
p=Geod(ellps='WGS84')
maxdist=0
for k in range(len(sta_lonlat)):
    lonsta=sta_lonlat[k,0]*ones(len(fault))
    latsta=sta_lonlat[k,1]*ones(len(fault))
    az,baz,d=p.inv(lonsta,latsta,fault[:,1],fault[:,2])
    d=d/1000 #to km
    if d.max()>maxdist:
        maxdist=d.max()
        
print 'Maximum distance is '+str(maxdist)+' km'
    


# Save new fault file
savetxt(fout,fault,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e',header='#')


##Make srf
forward.mudpy2srf(fout,log)
