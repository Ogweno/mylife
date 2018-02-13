from mudpy.hfsims import radiation_pattern
from numpy import arange,zeros,deg2rad,ones,linspace
from matplotlib import pyplot as plt

strike=0
dip=89
take_off_angle=45
rake=180

azimuth=arange(0,360)
P=zeros(len(azimuth))
SV=zeros(len(azimuth))
SH=zeros(len(azimuth))

for k in range(len(azimuth)):
    P[k],SV[k],SH[k]=radiation_pattern(strike,dip,rake,azimuth[k],take_off_angle)

plt.figure()    
ax = plt.subplot(131, projection='polar')
ax.plot(deg2rad(azimuth), abs(P))
ax = plt.subplot(132, projection='polar')
ax.plot(deg2rad(azimuth), abs(SV))
ax = plt.subplot(133, projection='polar')
ax.plot(deg2rad(azimuth), abs(SH))



N=100
strike=0*ones(N)
dip=89*ones(N)
take_off_angle=45*ones(N)
rake=180*ones(N)
azimuth=linspace(0,360,N)

P,SV,SH=radiation_pattern(strike,dip,rake,azimuth,take_off_angle)

plt.figure()    
ax = plt.subplot(131, projection='polar')
ax.plot(deg2rad(azimuth), abs(P))
ax = plt.subplot(132, projection='polar')
ax.plot(deg2rad(azimuth), abs(SV))
ax = plt.subplot(133, projection='polar')
ax.plot(deg2rad(azimuth), abs(SH))

plt.show()
