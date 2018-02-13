from numpy import genfromtxt,r_,where,mean,log10,array
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
import nitime.algorithms as tsa
from matplotlib.ticker import MultipleLocator
from mtspec import mt_coherence
from numpy import arange,zeros
from scipy.interpolate import interp1d
from run_filt import RunningMedian as rm
from obspy import UTCDateTime

#st=read('/Users/dmelgar/tidegauge_noise/eureka/eure_2017.sac')
#vmin=-1.5
#vmax=0.5
#sealim=0.15
#specmax=2.5

st=read(u'/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_2017.sac')
stp=read(u'/Users/dmelgar/tidegauge_noise/data/cres/sac/cres_pres_2017.mseed')

Tw=10
Ntapers=14

dt=int(round(st[0].stats.delta))



#trim
st[0].trim(starttime=UTCDateTime('20170301T00:00:00'))
stp.trim(starttime=UTCDateTime('20170301T00:00:00'))

#demean
st[0].data=st[0].data-mean(st[0].data)
stp[0].data=stp[0].data-mean(stp[0].data)

out = mt_coherence(dt, st[0].data, stp[0].data, Tw, Ntapers, st[0].stats.npts, 0.95, freq=True,cohe=True, iadapt=1)


plt.plot(1/out['freq']/60, out['cohe'])
plt.xlabel("Frequency [Hz]")
plt.ylabel("Coherency")
plt.xlim([0,120])
plt.show()