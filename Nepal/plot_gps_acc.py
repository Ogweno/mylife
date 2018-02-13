from obspy import read
from numpy import r_,diff
from matplotlib import pyplot as plt
import matplotlib

fcorner=0.5

n1=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXN.sac')
e1=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXE.sac')
u1=read(u'/Users/dmelgar/Nepal2015/GPS/cut/KKN4.LXZ.sac')

n2=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXN.sac')
e2=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXE.sac')
u2=read(u'/Users/dmelgar/Nepal2015/GPS/cut/NAST.LXZ.sac')

n3=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.n')
e3=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.e')
u3=read(u'/Users/dmelgar/Nepal2015/strong_motion/KATNP.acc.u')

matplotlib.rcParams.update({'font.size': 14})

fig, axarr = plt.subplots(3, 3)  
axn=axarr[0,0]
axe=axarr[1,0]
axu=axarr[2,0]

axn.plot(n1[0].times()-10+12,n1[0].data,'k')
axn.set_xlim([0,60])
axn.set_xticklabels([])
axn.grid()
axn.set_yticks([-2,0])
axn.set_ylabel('North')
axn.set_xticks([0,20,40,60])
axe.plot(e1[0].times()-10+12,e1[0].data,'k')
axe.set_xlim([0,60])
axe.set_xticklabels([])
axe.set_yticks([0,0.5])
axe.grid()
axe.set_ylabel('East')
axe.set_xticks([0,20,40,60])
axu.plot(u1[0].times()-10+12,u1[0].data,'k')
axu.set_xlim([0,60])
axu.set_ylim([-0.2,1.7])
axu.set_yticks([0,1.5])
axu.set_xticks([0,20,40,60])
axu.grid()
axu.set_xticks([0,20,40,60])
axu.set_ylabel('Up')
axn.set_title('KKN4 (m)',fontsize=14)


axn=axarr[0,1]
axe=axarr[1,1]
axu=axarr[2,1]

axn.plot(n2[0].times()-10+14.4,n2[0].data,'k')
axn.set_xlim([0,60])
axn.set_xticklabels([])
axn.grid()
axn.set_yticks([-1,0])
axn.set_xticks([0,20,40,60])
axe.plot(e2[0].times()-10+14.4,e2[0].data,'k')
axe.set_xlim([0,60])
axe.set_xticklabels([])
axe.set_yticks([0,0.5])
axe.grid()
axe.set_xticks([0,20,40,60])
axu.plot(u2[0].times()-10+14.4,u2[0].data,'k')
axu.set_xlim([0,60])
axu.set_ylim([-0.2,1.2])
axu.set_yticks([0,1])
axu.set_xticks([0,20,40,60])
axu.grid()
axu.set_xlabel('Seconds')
axu.set_xticks([0,20,40,60])
axn.set_title('NAST (m)',fontsize=14)

axn=axarr[0,2]
axe=axarr[1,2]
axu=axarr[2,2]

axn.plot(n3[0].times()-10+13.45,n3[0].data,'k')
axn.set_xlim([0,90])
axn.set_xticklabels([])
axn.grid()
axn.set_yticks([-1,1])
axn.set_xticks([0,20,40,60])
axe.plot(e3[0].times()-10+13.45,e3[0].data,'k')
axe.set_xlim([0,90])
axe.set_xticklabels([])
axe.set_yticks([-1,1])
axe.grid()
axe.set_xticks([0,20,40,60])
axu.plot(u3[0].times()-10+13.45,u3[0].data,'k')
axu.set_xlim([0,90])
axu.set_yticks([-1,1])
axu.set_xticks([0,30,60,90])
axu.grid()
axu.set_xticks([0,20,40,60])
axn.set_title(r'KATNP (m/s$^2$)',fontsize=14)

plt.subplots_adjust(left=0.2, bottom=0.15, right=0.8, top=0.85, wspace=0.25, hspace=0)


plt.show()


