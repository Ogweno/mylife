from numpy import genfromtxt,r_,where
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime

#Read actual data
crot_short=read(u'/Users/dmelgar/Lefkada2015/tsunami/sac/crot.sac')
time_epi=UTCDateTime('2015-11-17T07:10:07')
t1=crot_short[0].stats.starttime-time_epi

#read simulation
tg=genfromtxt('/Users/dmelgar/Tsunamis/Lefkada_regional_noadvection/_output/fort.gauge')
i=where(tg[:,0]==1003)[0]
t=tg[i,2]
crot=tg[i,6]
t=r_[-3600,t]
crot=r_[0,crot]

plt.figure(figsize=(6,3))

ax=plt.subplot(211)
plt.plot((crot_short[0].times()+t1)/3600,crot_short[0].data,c='#606060')
plt.plot([0,0],[-0.13,0.13],'--',c='k',lw=2)
ax.set_xlim([-1,4])
ax.set_ylim([-0.13,0.13])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','-0.1','','0','','0.1',''])
ax.set_ylabel('Tsunami (m)')

ax.annotate('CR08',xy=(-0.8,0.07),fontsize=14)
ax.annotate('obs',xy=(-0.8,-0.1),fontsize=14)

ax=plt.subplot(212)
plt.plot((t)/3600,crot,'r')
ax.set_xlim([-1,4])
ax.set_xlabel('Hours after OT')
ax.set_ylabel('Tsunami (m)')
ax.set_ylim([-0.003,0.003])
plt.plot([0,0],[-0.02,0.02],'--',c='k',lw=2)
#ax.yaxis.set_ticklabels(['-0.02','','-0.01','','0','','0.01','','0.02',''])
ax.yaxis.set_ticklabels(['','-0.002','','0','','0.002','',''])
ax.annotate('synth',xy=(-0.8,-0.0023),fontsize=14)

plt.subplots_adjust(hspace=0.03,bottom=0.2,right=0.98,left=0.15)

plt.show()