from obspy import read
from matplotlib import pyplot as plt
from numpy import arange,genfromtxt
from butter import bandpass

st1=read('/Users/dmelgar/Slip_inv/Napa_seis/structure/gil7.mod_15/30.syn.r')
st1[0].data=st1[0].data*1000
stf1=genfromtxt('/Users/dmelgar/Slip_inv/Napa_seis/structure/gil7.mod_15/floatArray.txt')
tstf1=arange(0,len(stf1)*0.2-0.2,0.2)

st1f=st1.copy()
st1f[0].data=bandpass(st1f[0].data,[1./20,1.0],0.2,2)

st2=read('/Users/dmelgar/Slip_inv/Napa_seis/structure/gil7hi.mod_15/30.syn.r')
st2[0].data=st2[0].data*1000
stf2=genfromtxt('/Users/dmelgar/Slip_inv/Napa_seis/structure/gil7hi.mod_15/floatArray.txt')
tstf2=arange(0,len(stf2)*0.02,0.02)

st2f=st2.copy()
st2f[0].data=bandpass(st2f[0].data,[1./20,1.0],0.02,2)


plt.figure()
plt.subplot(231)
plt.plot(st1[0].times()-st1[0].times()[49],st1[0].data)
plt.title('5Hz with what we thought was 0.1s rise time STF')
plt.xlim([-1,30])
plt.ylabel('d(mm)')
plt.xlabel('Time (s)')

plt.subplot(232)
plt.plot(st1f[0].times()-st1f[0].times()[49],st1f[0].data)
plt.title('Bandpass fitlered between 1/20 to 1 Hz')
plt.xlim([-1,30])
plt.xlabel('Time (s)')

plt.subplot(233)
plt.plot(tstf1,stf1)
plt.title('Actual STF used @ 5Hz')
plt.ylabel('Slip rate (m/s)')
plt.xlabel('Time (s)')

plt.subplot(234)
plt.plot(st2[0].times()-st2[0].times()[49],st2[0].data)
plt.title('50Hz with actually 0.1s rise time STF')
plt.xlim([-1,30])
plt.ylabel('d(mm)')
plt.xlabel('Time (s)')

plt.subplot(235)
plt.plot(st2f[0].times()-st2[0].times()[49],st2f[0].data)
plt.title('Bandpass fitlered between 1/20 to 1 Hz')
plt.xlim([-1,30])
plt.xlabel('Time (s)')

plt.subplot(236)
plt.plot(tstf2,stf2)
plt.title('Actual STF used @ 50Hz')
plt.ylabel('Slip rate (m/s)')
plt.xlabel('Time (s)')

plt.figure()
plt.plot(st1f[0].times()-st1[0].times()[49],st1f[0].data)
plt.plot(st2f[0].times()-st2[0].times()[49],st2f[0].data)
plt.title('Bandpass fitlered between 1/20 to 1 Hz')
plt.xlim([-1,30])
plt.ylabel('d(m)')
plt.xlabel('Time (s)')

plt.show()