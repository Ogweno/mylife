from numpy import genfromtxt,zeros,c_,mean,arange,where
from obspy.core import UTCDateTime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mudpy.forward import lowpass
from datetime import timedelta
from obspy import Stream,Trace
from obspy import read


st1=read('/Users/dmelgar/Chiapas2017/tsunami/long/ptan.sac')
st2=read('/Users/dmelgar/Chiapas2017/tsunami/long/huat.sac')
st3=read('/Users/dmelgar/Chiapas2017/tsunami/long/sali.sac')
st4=read('/Users/dmelgar/Chiapas2017/tsunami/long/pchi.sac')


xl=[-24,72]
tminus=24

plt.figure(figsize=(12,4))


plt.subplot(411)
plt.plot(st1[0].times()/3600-tminus,st1[0].data)
plt.ylabel('PTAN')
plt.xlim(xl)


plt.subplot(412)
plt.plot(st2[0].times()/3600-tminus,st2[0].data)
plt.ylabel('HUAT')
plt.xlim(xl)

plt.subplot(413)
plt.plot(st3[0].times()/3600-tminus,st3[0].data)
plt.ylabel('SALI')
plt.xlim(xl)

plt.subplot(414)
plt.plot(st4[0].times()/3600-tminus,st4[0].data)
plt.ylabel('PMAD')
plt.xlim(xl)
plt.xlabel('Hours since OT')

plt.show()