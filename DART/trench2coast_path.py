from numpy import genfromtxt,where,ones,zeros,argmin,savetxt,linspace
from matplotlib import pyplot as plt
from pyproj import Geod

coast=genfromtxt('/Users/dmelgar/DART_analysis/coasts.txt')
trench=genfromtxt(u'/Users/dmelgar/DART_analysis/all_trenches.txt')
fout_coast=u'/Users/dmelgar/DART_analysis/nearest_coast.txt'
fout_profiles=u'/Users/dmelgar/DART_analysis/all_profiles.txt'

#Get rid of negative lons
#i=where(coast[:,0]<0)[0]
#coast[i,0]=360+coast[i,0]
#
#i=where(trench[:,0]<0)[0]
#trench[i,0]=360+trench[i,0]

coast_out=zeros((len(trench),2))
p=Geod(ellps='WGS84')
for k in range(len(trench)):
    print k
    az,baz,dist=p.inv(coast[:,0],coast[:,1],ones(len(coast))*trench[k,0],ones(len(coast))*trench[k,1])
    i=argmin(dist)
    coast_out[k,:]=coast[i,:]

savetxt(fout_coast,coast_out,fmt='%.4f')



#Get profiles
Npts=20
profiles=zeros((len(trench)*2,Npts))
for k in range(len(profiles)/2):
    print k
    x=linspace(coast_out[k,0],trench[k,0],Npts)
    y=linspace(coast_out[k,1],trench[k,1],Npts)
    profiles[2*k,:]=x
    profiles[2*k+1,:]=y

savetxt(fout_profiles,profiles,fmt='%.4f')


plt.figure()
plt.scatter(coast[:,0],coast[:,1],s=3)
plt.scatter(trench[:,0],trench[:,1],s=3)