from obspy import read
from numpy import genfromtxt,where,r_
from matplotlib import pyplot as plt
import matplotlib as mpl
label_size = 14
mpl.rcParams['ytick.labelsize'] = label_size 


plt.figure(figsize=(4,12))


g=genfromtxt(u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/gauge01245.txt')
i=where(g[:,0]>=2)[0]
data=read(u'/Users/dmelgar/NewZealand2016/tsunami/sac/cast.sac')
xl=[0,240]
yl=[-0.25,0.25]
ax=plt.subplot(811)
ax.plot([-10,-10],[0,0],c='#1E90FF')
ax.plot([-10,-10],[0,0],'#FF6347')
ax.legend(['obs','syn'],ncol=2,fontsize=14,loc=1,bbox_to_anchor=(0.9, 1.7),frameon=False)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(data[0].times()/60,data[0].data,'#1E90FF')
ax.plot([0,240],[0,0],lw=0.5)
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.annotate(s='cpit',xy=(2,0.1),fontsize=16)
ax.set_ylabel('m')
ax=plt.subplot(812)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(r_[0,g[i,1]/60],r_[0,g[i,5]],'#FF6347')
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.set_ylabel('m')
ax.tick_params(labelbottom='off')




g=genfromtxt(u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/gauge01224.txt')
i=where(g[:,0]>=2)[0]
data=read(u'/Users/dmelgar/NewZealand2016/tsunami/sac/well.sac')
xl=[0,240]
yl=[-0.4,0.4]
ax=plt.subplot(813)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(data[0].times()/60,data[0].data,'#1E90FF')
ax.plot([0,240],[0,0],lw=0.5)
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.annotate(s='wlgt',xy=(2,0.1),fontsize=16)
ax.set_ylabel('m')
ax=plt.subplot(814)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(r_[0,g[i,1]/60],r_[0,g[i,5]],'#FF6347')
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.set_ylabel('m')



g=genfromtxt(u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/gauge01149.txt')
data=genfromtxt(u'/Users/dmelgar/code/GMT/NewZealand/kaik.txt')
xl=[0,240]
yl=[-3,1.5]
ax=plt.subplot(815)

ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot([0,240],[-0.89,-0.89],'--',lw=0.5,c='k')
ax.plot(data[:,0],data[:,1],'#1E90FF')
ax.plot([0,240],[0,0],lw=0.5)
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.annotate(s='kait',xy=(2,0.2),fontsize=16)
ax.set_ylabel('m')
ax=plt.subplot(816)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot([0,240],[-0.89,-0.89],'--',lw=0.5,c='k')
ax.plot(g[3:,1]/60,g[3:,5]-0.89,'#FF6347')
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.set_ylabel('m')


g=genfromtxt(u'/Users/dmelgar/Tsunamis/Kaikoura_gfast/_output/gauge01093.txt')
i=where((g[:,0]>=3) & (g[:,5]<1.0))[0]
data=read(u'/Users/dmelgar/NewZealand2016/tsunami/sac/sumt.sac')
xl=[0,240]
yl=[-0.7,0.7]
ax=plt.subplot(817)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(data[0].times()/60,data[0].data/1e5,'#1E90FF')
ax.plot([0,240],[0,0],lw=0.5)
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.tick_params(labelbottom='off')
ax.annotate(s='sumt',xy=(2,0.25),fontsize=16)
ax.set_ylabel('m')
ax=plt.subplot(818)
ax.plot([0,240],[0,0],lw=0.5,c='k')
ax.plot(g[i,1]/60,g[i,5],'#FF6347')
ax.set_ylim(yl)
ax.set_xlim(xl)
ax.set_xlabel('Minutes after OT',fontsize=14)
ax.set_ylabel('m')

plt.subplots_adjust(left=0.25,right=0.98,hspace=0.3)

plt.show()

