from obspy import read
from datetime import timedelta
from matplotlib import pyplot as plt
from numpy import r_,ones,diff
from scipy.signal import deconvolve

N=12
tcut=10
tend=60

kkn4=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXE.sac')
nast=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXE.sac')
kkn4[0].data=r_[ones(N)*kkn4[0].data[0],kkn4[0].data]

kkn4.trim(starttime=kkn4[0].stats.starttime+timedelta(seconds=tcut),endtime=kkn4[0].stats.starttime+timedelta(seconds=tend))
nast.trim(starttime=nast[0].stats.starttime+timedelta(seconds=tcut),endtime=nast[0].stats.starttime+timedelta(seconds=tend))

kkn4[0].data=r_[1e-5,diff(kkn4[0].data)/0.2]
nast[0].data=r_[1e-5,diff(nast[0].data)/0.2]

basin,remainder=deconvolve(kkn4[0].data,nast[0].data)

plt.figure()
plt.plot(kkn4[0].times(),kkn4[0].data,nast[0].times(),nast[0].data)
plt.show()