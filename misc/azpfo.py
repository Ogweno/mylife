from numpy import genfromtxt,ones
from pyproj import Geod
from matplotlib import pyplot as plt

times=genfromtxt('/Users/dmelgar/Downloads/AZ_PFO_catalogue.txt',usecols=[0],dtype='S')
data=genfromtxt('/Users/dmelgar/Downloads/AZ_PFO_catalogue.txt',usecols=[1,2,3,4])
station_lon=ones(len(data))*(-116.459400)
station_lat=ones(len(data))*33.611700
g=Geod(ellps='WGS84')
az,baz,dist=g.inv(station_lon,station_lat,data[:,1],data[:,0])
dist=dist/1000.
f=open('/Users/dmelgar/Downloads/AZ_PFO_catalogue_500km.txt','w')
f.write('# Event time                     lon             lat     depth(km)   Mag    Dist(km)\n')
for k in range(len(data)):
    line='%s\t%12.6f\t%12.6f\t%.2f\t%.2f\t%.2f\n' % (times[k],data[k,1],data[k,0],data[k,2],data[k,3],dist[k])
    f.write(line)
f.close()

plt.subplot(121)
plt.scatter(data[:,1],data[:,0])
plt.scatter(-116.459400,33.611700,marker='*',s=450,c='r')
plt.xlabel('Lon')
plt.ylabel('Lat')

plt.subplot(122)
plt.plot(data[:,3])
plt.ylabel('Magnitude')
plt.xlabel('Event No.')

