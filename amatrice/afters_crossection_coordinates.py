from pyproj import Geod
from numpy import genfromtxt,r_,zeros
from matplotlib import pyplot as plt

#origin=[13.1031,42.4078]
origin=[13.04139759,42.50569016]
strike=155+180.
dx=5e3 #in m
prof_length=45e3 #in m
nprofiles=10


af1=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records.txt',usecols=[0,1,2,3,4,5])
af2=genfromtxt('/Users/dmelgar/Amatrice2016/afters/basile/hypodd_ct_records2.txt',usecols=[0,1,2,3,4,5])
afters=r_[af1,af2]

p=Geod(ellps='WGS84')

        

profiles=zeros((nprofiles,4))
for k in range(nprofiles):
    #Get start of profile
    profiles[k,0],profiles[k,1],baz=p.fwd(origin[0],origin[1],strike,dx*(k+1))
    #Get end
    profiles[k,2],profiles[k,3],baz=p.fwd(profiles[k,0],profiles[k,1],strike+90,prof_length)
    
    
#COORDINATES OF 10,20 30 AND 40KM LINES
line10=zeros((2,2))
line10[0,0],line10[0,1],baz=p.fwd(profiles[0,0],profiles[0,1],strike+90,10e3)
line10[1,0],line10[1,1],baz=p.fwd(profiles[-1,0],profiles[-1,1],strike+90,10e3)
line20=zeros((2,2))
line20[0,0],line20[0,1],baz=p.fwd(profiles[0,0],profiles[0,1],strike+90,20e3)
line20[1,0],line20[1,1],baz=p.fwd(profiles[-1,0],profiles[-1,1],strike+90,20e3)
line30=zeros((2,2))
line30[0,0],line30[0,1],baz=p.fwd(profiles[0,0],profiles[0,1],strike+90,30e3)
line30[1,0],line30[1,1],baz=p.fwd(profiles[-1,0],profiles[-1,1],strike+90,30e3)
line40=zeros((2,2))
line40[0,0],line40[0,1],baz=p.fwd(profiles[0,0],profiles[0,1],strike+90,40e3)
line40[1,0],line40[1,1],baz=p.fwd(profiles[-1,0],profiles[-1,1],strike+90,40e3)

    
    
plt.figure()
plt.scatter(afters[:,3],afters[:,2],lw=0,c='b')
for k in range(nprofiles):
    plt.plot(profiles[k,[0,2]],profiles[k,[1,3]],'k')
plt.plot(line10[:,0],line10[:,1])
plt.plot(line20[:,0],line20[:,1])
plt.plot(line30[:,0],line30[:,1])
plt.plot(line40[:,0],line40[:,1])
plt.axis('equal')
plt.show()