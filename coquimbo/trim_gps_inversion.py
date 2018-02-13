from matplotlib import pyplot as plt
from obspy import read
from datetime import timedelta

sta='CNBA'
path=u'/Users/dmelgar/Coquimbo2015/GPS/trim/120s/'
outpath='/Users/dmelgar/Slip_inv/Coquimbo_4s/data/waveforms/'
tcut=timedelta(seconds=90)

n=read(path+sta+'.LXN.sac')
e=read(path+sta+'.LXE.sac')
u=read(path+sta+'.LXZ.sac')

plt.close("all")
plt.figure()

plt.subplot(321)
plt.plot(n[0].times(),n[0].data)
plt.xlim([0,130])
plt.subplot(323)
plt.plot(e[0].times(),e[0].data)
plt.xlim([0,130])
plt.subplot(325)
plt.plot(u[0].times(),u[0].data)
plt.xlim([0,130])

n.trim(endtime=n[0].stats.starttime+tcut)
e.trim(endtime=e[0].stats.starttime+tcut)
u.trim(endtime=u[0].stats.starttime+tcut)

plt.subplot(322)
plt.plot(n[0].times(),n[0].data)
plt.xlim([0,130])
plt.subplot(324)
plt.plot(e[0].times(),e[0].data)
plt.xlim([0,130])
plt.subplot(326)
plt.plot(u[0].times(),u[0].data)
plt.xlim([0,130])
plt.show()

n.write(outpath+sta+'.disp.n',format='SAC')
e.write(outpath+sta+'.disp.e',format='SAC')
u.write(outpath+sta+'.disp.u',format='SAC')