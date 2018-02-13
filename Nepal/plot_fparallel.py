from obspy import read,Stream,Trace
from matplotlib import pyplot as plt
from numpy import diff
r=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXR.sac')
t=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXT.sac')
z=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXZ.sac')
vz=Stream(Trace())
vz=z.copy()
vz[0].data=diff(vz[0].data)/0.2

plt.figure()
plt.plot(r[0].times(),r[0].data,t[0].times(),t[0].data,z[0].times(),z[0].data,vz[0].times(),vz[0].data)
plt.annotate('Fault normal displ.(m)',xy=(38,-2.2))
plt.annotate('Fault parallel displ.(m)',xy=(38,-0.6))
plt.annotate('Vert. displacement(m)',xy=(38,1.0))
plt.annotate('Vert. velocity(m/s)',xy=(38,0.12))
plt.plot([20,20],[-10,10],'k--')
plt.plot([28,28],[-10,10],'k--')
plt.ylim([-2.5,2])
plt.xlim([0,60])
plt.xlabel('Seconds',fontsize=14)
plt.show()