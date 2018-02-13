from numpy import genfromtxt,r_,where
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




#Read actual data
crot_short=read(u'/Users/dmelgar/Lefkada2015/tsunami/sac/crot.sac')
time_epi=UTCDateTime('2015-11-17T07:10:07')
t1=crot_short[0].stats.starttime-time_epi

#read simulation
tg=read('/Users/dmelgar/Proposals/NEHRP_2016/indonesia/meulaboh_2012.sac')

plt.figure(figsize=(6,6))

ax=plt.subplot(211)
plt.plot((crot_short[0].times()+t1)/3600,crot_short[0].data,c='#1E90FF',lw=2)
ax.set_xlim([0,4.5])
ax.set_ylim([-0.13,0.13])
ax.set_ylabel('Tsunami (m)')
xmajorLocator = MultipleLocator(1)
xminorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.1)
yminorLocator = MultipleLocator(0.02)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.tick_params(which='major',length=7,width=1)
ax.tick_params(which='minor',length=4,width=1)
plt.xlabel('Hours after OT')


ax=plt.subplot(212)
plt.plot((tg[0].times())/3600,tg[0].data,c='#1E90FF',lw=2)
ax.set_xlim([0,4.5])
ax.set_ylim([-0.55,1.1])
ax.set_ylabel('Tsunami (m)')
xmajorLocator = MultipleLocator(1)
xminorLocator = MultipleLocator(0.1)
ymajorLocator = MultipleLocator(0.5)
yminorLocator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.tick_params(which='major',length=7,width=1)
ax.tick_params(which='minor',length=4,width=1)
plt.xlabel('Hours after OT')

#ax.annotate('CR08',xy=(-0.8,0.07),fontsize=14)
#ax.annotate('obs',xy=(-0.8,-0.1),fontsize=14)


plt.subplots_adjust(hspace=0.4,bottom=0.2,right=0.98,left=0.15)

plt.show()