from matplotlib import pyplot as plt
from numpy import genfromtxt
from scipy.integrate import cumtrapz

v1=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.nwhp.vel.bbp')
v2=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.vasq.vel.bbp')
v3=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.balp.vel.bbp')

a1=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.nwhp.acc.bbp')
a2=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.vasq.acc.bbp')
a3=genfromtxt('/Users/dmelgar/code/BBP/bbp/bbp_data/outdata/7301769/7301769.balp.acc.bbp')

#Integrate
v1[:,1]=cumtrapz(v1[:,1],v1[:,0],initial=0)
v1[:,2]=cumtrapz(v1[:,2],v1[:,0],initial=0)
v1[:,3]=cumtrapz(v1[:,3],v1[:,0],initial=0)
v2[:,1]=cumtrapz(v2[:,1],v2[:,0],initial=0)
v2[:,2]=cumtrapz(v2[:,2],v2[:,0],initial=0)
v2[:,3]=cumtrapz(v2[:,3],v2[:,0],initial=0)
v3[:,1]=cumtrapz(v3[:,1],v3[:,0],initial=0)
v3[:,2]=cumtrapz(v3[:,2],v3[:,0],initial=0)
v3[:,3]=cumtrapz(v3[:,3],v3[:,0],initial=0)


plt.figure()
yl=[-5,5]
xl=[0,40]

plt.subplot(321)
plt.title('Vertical acceleration (cm/s/s)')
plt.plot(a1[:,0],a1[:,3])
plt.ylim(yl)
plt.xlim(xl)
plt.ylabel('NHWP')

plt.subplot(323)
plt.plot(a2[:,0],a2[:,3])
plt.ylim(yl)
plt.xlim(xl)
plt.ylabel('VASQ')

plt.subplot(325)
plt.plot(a3[:,0],a3[:,3])
plt.ylim(yl)
plt.xlim(xl)
plt.ylabel('BALP')
plt.xlabel('Seconds after OT')


yl=[-1,5]
plt.subplot(322)
plt.title('Vertical displacement (cm)')
plt.plot(v1[:,0],v1[:,3])
plt.ylim(yl)
plt.xlim(xl)

yl=[-3,1]
plt.subplot(324)
plt.plot(v2[:,0],v2[:,3])
plt.ylim(yl)
plt.xlim(xl)

yl=[-0.1,0.1]
plt.subplot(326)
plt.plot(v3[:,0],v3[:,3])
plt.ylim(yl)
plt.xlim(xl)
plt.xlabel('Seconds after OT')




#plt.figure()
#
#plt.subplot(311)
#plt.plot(v1[:,0],v1[:,1])
#plt.ylabel('North (cm)')
#plt.title('NWHP displacement')
#plt.xlim([0,60])
#
#plt.subplot(312)
#plt.plot(v1[:,0],v1[:,2])
#plt.ylabel('East (cm)')
#plt.xlim([0,60])
#
#plt.subplot(313)
#plt.plot(v1[:,0],v1[:,3])
#plt.ylabel('Up (cm)')
#plt.xlim([0,60])
#
#plt.figure()
#
#plt.subplot(311)
#plt.plot(v2[:,0],v2[:,1])
#plt.ylabel('North (cm)')
#plt.title('BALP displacement')
#plt.xlim([0,60])
#
#plt.subplot(312)
#plt.plot(v2[:,0],v2[:,2])
#plt.ylabel('East (cm)')
#plt.xlim([0,60])
#
#plt.subplot(313)
#plt.plot(v2[:,0],v2[:,3])
#plt.ylabel('Up (cm)')
#plt.xlim([0,60])
#
#plt.figure()
#
#plt.subplot(311)
#plt.plot(v3[:,0],v3[:,1])
#plt.ylabel('North (cm)')
#plt.title('VASQ displacement')
#plt.xlim([0,60])
#
#plt.subplot(312)
#plt.plot(v3[:,0],v3[:,2])
#plt.ylabel('East (cm)')
#plt.xlim([0,60])
#
#plt.subplot(313)
#plt.plot(v3[:,0],v3[:,3])
#plt.ylabel('Up (cm)')
#plt.xlim([0,60])

plt.show()
