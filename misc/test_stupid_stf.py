from obspy import read
from scipy.signal import convolve
from mudpy.forward import build_source_time_function
from matplotlib import pyplot as plt

#synth='/Users/dmelgar/FakeQuakes/Hayward/GFs/misc/gil7.mod_5.5010.sub0017_5Hz/P224.subfault0017.SS.disp.e'
synth='/Users/dmelgar/FakeQuakes/Hayward/GFs/misc/gil7.mod_5.5010.sub0017_1Hz/P224.subfault0017.SS.disp.e'
#synth='/Users/dmelgar/FakeQuakes/Hayward/GFs/misc/gil7.mod_5.5010.sub0017_5Hz/15.270.grn.4'
#synth='/Users/dmelgar/FakeQuakes/Hayward/GFs/misc/gil7.mod_5.5010.sub0017_1Hz/15.270.grn.4'

st=read(synth)
sti=read(synth)

dt=st[0].stats.delta
rise_time=1.0
NFFT=st[0].stats.npts
total_time=NFFT*dt

#First build stupid STF and convolve WITHOUT interpolation
t1,Mdot1=build_source_time_function(rise_time,dt,total_time,stf_type='dreger',scale=True)
after_conv=convolve(st[0].data,Mdot1)[0:NFFT]


#Simple interpolationto 10Hz
sti[0].interpolate(10,method='linear')
dt=sti[0].stats.delta
NFFT=sti[0].stats.npts
total_time=dt*NFFT
#Now build an STF at 10Hz
t2,Mdot2=build_source_time_function(rise_time,dt,total_time,stf_type='dreger',scale=True)
#Convolution WITH interpoalted synthetic
after_conv_interp=convolve(sti[0].data,Mdot2)[0:NFFT]


plt.figure()
plt.subplot(211)
plt.plot(t1,Mdot1,'k')
plt.plot(t2,Mdot2,'r')
plt.xlim([0,2])
plt.xlabel('Seconds')
plt.legend(['Original','10Hz'])

plt.subplot(212)
plt.plot(st[0].times(),after_conv,'k')
plt.plot(sti[0].times(),after_conv_interp,'r')
plt.xlim([40,70])
plt.xlabel('Seconds')
plt.legend(['Original','10Hz'])

plt.show()