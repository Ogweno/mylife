from numpy import genfromtxt,log10,ones,savetxt
from pyproj import Geod

vr=2.6
hypo=[-101.819,17.598,20.0] 
fout='/Users/dmelgar/Michoacan1985/slip/mich_after_kinematic.rupt'
f=genfromtxt('/Users/dmelgar/Michoacan1985/slip/mich_after.rupt')

#Moment
mu=45e9
M0=f[:,9]*f[:,10]*f[:,11]*mu
M0=M0.sum()
Mw=(2./3)*(log10(M0)-9.1)

#Average rise time
Ta=4.308e-7*M0**(1./3)

trise=f[:,9]**0.5
k=Ta/trise.mean()
trise=trise*k

#rupture onset
p=Geod(ellps='WGS84')
az,baz,Hdist=p.inv(f[:,1],f[:,2],ones(len(f))*hypo[0],ones(len(f))*hypo[1])
Hdist=Hdist/1000.
Vdist=abs(f[:,3]-hypo[2])
dist=(Hdist**2+Vdist**2)**0.5
tonset=dist/vr

#finalize
f[:,7]=trise
f[:,12]=tonset

#1	-101.785924	 17.346399	  7.794857	 300.00	  14.00	 0.5000	 1.0000	 0.0000	 1.2500	15000.0000	13900.0000	 1.0000

savetxt(fout,f,fmt='%d\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.3f',header='#')