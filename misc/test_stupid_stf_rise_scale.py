from obspy import read
from scipy.signal import convolve
from mudpy.forward import build_source_time_function
from matplotlib import pyplot as plt

synth='/Users/dmelgar/FakeQuakes/Hayward/GFs/misc/gil7.mod_5.5010.sub0017_5Hz/P224.subfault0017.SS.disp.e'
st=read(synth)

rise_scale1=4
rise_scale2=8
dt=st[0].stats.delta
rise_time=1.5
NFFT=1024
total_time=NFFT*dt

#First build stupid STF and convolve with standard rise_scale=4
t1,Mdot1=build_source_time_function(rise_time,dt,total_time,stf_type='dreger',scale=True,rise_scale=rise_scale1)
after_conv1=convolve(st[0].data,Mdot1)[0:NFFT]

#First build stupid STF and convolve with  rise_scale=6
t2,Mdot2=build_source_time_function(rise_time,dt,total_time,stf_type='dreger',scale=True,rise_scale=rise_scale2)
after_conv2=convolve(st[0].data,Mdot2)[0:NFFT]


plt.figure()
plt.subplot(211)
plt.plot(t1,Mdot1,'k')
plt.plot(t2,Mdot2,'r')
plt.xlim([0,2])
plt.xlabel('Seconds')
plt.legend(['rise_scale='+str(rise_scale1),'rise_scale='+str(rise_scale2)])

plt.subplot(212)
plt.plot(st[0].times(),after_conv1,'k')
plt.plot(st[0].times(),after_conv2,'r')
plt.xlim([5,25])
plt.xlabel('Seconds')
plt.legend(['rise_scale='+str(rise_scale1),'rise_scale='+str(rise_scale2)])

plt.show()