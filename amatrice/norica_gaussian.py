from numpy import genfromtxt,argmax,argmin,ones,tile,squeeze,sin,deg2rad,exp,zeros,r_,arange,loadtxt,log10,c_,savetxt
from pyproj import Geod
from matplotlib import pyplot as plt
from mudpy import forward

f=genfromtxt(u'/Users/dmelgar/Amatrice2016/Norcia1979/norcia_east.fault')
A=1
sigmax=20.0 #Striek direction
sigmay=33.0 #Dip direction
strike=320
dip=70

#Get along strike and along dip distances from top left of model
i=argmin(f[:,2])
lat0=f[i,2]
lon0=f[i,1]

P=Geod(ellps='WGS84')

#Get distances of first row to origin, then tile for all rows
az,baz,dist=P.inv(ones(10)*lon0,ones(10)*lat0,f[0:10,1],f[0:10,2])
dist=dist/1000.
Dstrike=squeeze(tile(dist,(1,12)))

#Downdip
z0=f[0,3]
deltaz=abs(z0-f[:,3])
Ddip=deltaz/sin(deg2rad(dip))

#means
mux=Dstrike.max()/2
muy=Ddip.max()/2

#Make slip
slip=exp(-((Dstrike-mux)**2/sigmax+(Ddip-muy)**2/sigmay))

#Get moment of the thing                          GLOBALS                             ########
home='/Users/dmelgar/Slip_inv/'
project_name='Amatrice_3Dfitsgeol_final1'
fault_name='norcia_east.fault'
model_name='aci.mod' 
mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
mu=zeros(len(slip))
for k in range(len(slip)):
    mu[k]=forward.get_mu(mod,f[k,3])

A=f[:,-1]**2
#Compute moments
try:
    M0=mu*A*slip[:,0]
except:
    M0=mu*A*slip
#Total up and copute magnitude
M0=M0.sum()
Mw=(2./3)*(log10(M0)-9.1)

print 'Moment before scaling is '+str(M0)

#Scale
scale=M0/(7.32*10**17)
slip=slip/scale

try:
    M0=mu*A*slip[:,0]
except:
    M0=mu*A*slip
#Total up and copute magnitude
M0=M0.sum()
Mw=(2./3)*(log10(M0)-9.1)

print 'Moment after scaling is '+str(M0)

#Save to file
out=c_[f[:,0:8],zeros(len(slip)),-slip,f[:,8:10],zeros(len(slip)),mu]
fout=u'/Users/dmelgar/Amatrice2016/Norcia1979/norcia_east.rupt'
fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
savetxt(fout,out,fmt=fmtout)

#To coulomb
forward.inv2coulomb(u'/Users/dmelgar/Amatrice2016/Norcia1979/norcia_east.rupt',[13.3,42.6],'/Users/dmelgar/Amatrice2016/Norcia1979/norcia_east.coul')

plt.figure()
plt.scatter(f[:,1],f[:,2],c=slip,s=80,lw=0)
plt.colorbar()
plt.show()