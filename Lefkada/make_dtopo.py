from numpy import genfromtxt,r_,c_,zeros,ones,savetxt,where,array,arange,where
from obspy.core.util.geodetics import gps2DistAzimuth

#f='/Users/dmelgar/Lefkada2015/scripts/up_coseis.xy'
#fout='/Users/dmelgar/Lefkada2015/scripts/up_coseis.dtopo'
f='/Users/dmelgar/Lefkada2015/scripts/total_up.xy'
fout='/Users/dmelgar/Lefkada2015/scripts/total_up.dtopo'

dt=1.0
vs=3300.
epi=array([20.6002,38.6655])

#load data
d=genfromtxt(f)
lonlat=d[:,0:2]  #lon and lat
d=d[:,2] #Actual displacememnt
#Compute distance from hypocenter to all points
dist=zeros(len(d))
for k in range(len(d)):
    if k%100==0:
        print str(k)+'/'+str(len(d))
    dist[k],az,baz=gps2DistAzimuth(epi[1],epi[0],lonlat[k,1],lonlat[k,0])
#Form full length dtopo vector
tmax=dist.max()/vs
tvec=arange(0,tmax+dt,dt)
t=array([])
#Form time vectopr
for k in range(len(tvec)):
    t=r_[t,ones(len(d))*tvec[k]]
#Form coordinate and dispalcememnt vector
dout=array([])
for k in range(len(tvec)):
    dmax=tvec[k]*vs
    print dmax
    i=where(dist<dmax)[0]
    mask=zeros(len(d))
    mask[i]=1
    print len(i)
    dout=r_[dout,d*mask]

#Form lonlat vector
L=lonlat
for k in range(len(tvec)-1):
    L=r_[L,lonlat]

#The total
dtopo=c_[t,L,dout]
savetxt(fout,dtopo,fmt='%i\t%.6f\t%.6f\t%.6e')