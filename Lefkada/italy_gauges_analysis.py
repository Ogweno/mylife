from numpy import genfromtxt,r_,where
from obspy import read
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime

#Read actual data
crot_short=read(u'/Users/dmelgar/Lefkada2015/tsunami/sac/crot.sac')
otra_short=read(u'/Users/dmelgar/Lefkada2015/tsunami/sac/otra.sac')
tara_short=read(u'/Users/dmelgar/Lefkada2015/tsunami/sac/tara.sac')
time_epi=UTCDateTime('2015-11-17T07:10:07')
t1=crot_short[0].stats.starttime-time_epi

#read simulation
tg=genfromtxt('/Users/dmelgar/Tsunamis/Lefkada_regional_noadvection/_output/fort.gauge')
i=where(tg[:,0]==1003)[0]
t=tg[i,2]
crot=tg[i,6]
t=r_[-3600,t]
crot=r_[0,crot]*100

i=where(tg[:,0]==1007)[0]
totra=tg[i,2]
otra=tg[i,6]
totra=r_[-3600,totra]
otra=r_[0,otra]*100

i=where(tg[:,0]==1006)[0]
ttara=tg[i,2]
tara=tg[i,6]
ttara=r_[-3600,ttara]
tara=r_[0,tara]*100
i=where(ttara<1800)[0]
tara[i]=0




plt.figure(figsize=(6,9))

ax=plt.subplot(611)
plt.plot((crot_short[0].times()+t1)/3600,crot_short[0].data*100,c='#606060')
plt.plot([0,0],[-13,13],'--',c='k',lw=2)
ax.set_xlim([-1,4])
ax.set_ylim([-13,13])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','-10.0','','0','','10.0',''])

ax.annotate('CR08',xy=(-0.8,7),fontsize=14)
ax.annotate('obs',xy=(-0.8,-10),fontsize=14)

ax=plt.subplot(612)
plt.plot((t)/3600,crot,'r')
ax.set_xlim([-1,4])
ax.set_ylim([-0.3,0.3])
plt.plot([0,0],[-0.4,0.4],'--',c='k',lw=2)
#ax.yaxis.set_ticklabels(['-0.02','','-0.01','','0','','0.01','','0.02',''])
ax.yaxis.set_ticklabels(['','-0.2','','0','','0.2','',''])
ax.annotate('synth',xy=(-0.8,-0.23),fontsize=14)


ax=plt.subplot(613)
plt.plot((tara_short[0].times()+t1)/3600,tara_short[0].data*100,c='#606060')
plt.plot([0,0],[-13,13],'--',c='k',lw=2)
ax.set_xlim([-1,4])
ax.set_ylim([-5,5])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','-5.0','','0','','5.0',''])

ax.annotate('TA18',xy=(-0.8,2.8),fontsize=14)
ax.annotate('obs',xy=(-0.8,-3.6),fontsize=14)

ax=plt.subplot(614)
plt.plot((ttara)/3600,tara,'r')
ax.set_xlim([-1,4])
ax.set_ylabel('Tsunami Amplitude (cm)')
ax.set_ylim([-0.3,0.3])
plt.plot([0,0],[-0.4,0.4],'--',c='k',lw=2)
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','-0.2','','0','','0.2','',''])
ax.annotate('synth',xy=(-0.8,-0.23),fontsize=14)

ax=plt.subplot(615)
plt.plot((otra_short[0].times()+t1)/3600,otra_short[0].data*100,c='#606060')
plt.plot([0,0],[-13,13],'--',c='k',lw=2)
ax.set_xlim([-1,4])
ax.set_ylim([-5,5])
ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels(['','-5.0','','0','','5.0',''])

ax.annotate('OT15',xy=(-0.8,2.8),fontsize=14)
ax.annotate('obs',xy=(-0.8,-3.6),fontsize=14)

ax=plt.subplot(616)
plt.plot((totra)/3600,otra,'r')
ax.set_xlim([-1,4])
ax.set_xlabel('Hours after OT')
ax.set_ylim([-0.3,0.3])
plt.plot([0,0],[-0.4,0.4],'--',c='k',lw=2)
#ax.yaxis.set_ticklabels(['-0.02','','-0.01','','0','','0.01','','0.02',''])
ax.yaxis.set_ticklabels(['','-0.2','','0','','0.2','',''])
ax.annotate('synth',xy=(-0.8,-0.23),fontsize=14)




plt.subplots_adjust(hspace=0.03,bottom=0.2,right=0.98,left=0.17)

plt.show()