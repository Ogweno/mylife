from numpy import genfromtxt,sqrt,where,array,savetxt
from matplotlib import pyplot as plt

g=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/static/unr_all.txt',skip_header=2)
sta=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/static/unr_all.txt',usecols=0,dtype='S',skip_header=2)
thresh=0.02

d=sqrt(g[:,3]**2+g[:,4]**2+g[:,5]**2)

i=where(d>thresh)[0]

#plt.quiver(g[i,1],g[i,2],g[i,3],g[i,4])
plt.quiver(g[:,1],g[:,2],g[:,3],g[:,4])
for k in range(len(sta)):
    plt.annotate(s=sta[k],xy=(g[k,1],g[k,2]))
plt.show()


#Make neu files
g=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/static/unr_filtered.txt')
sta=genfromtxt('/Users/dmelgar/Chiapas2017/GPS/static/unr_filtered.txt',usecols=0,dtype='S')

for k in range(len(sta)):
    neu=array([[g[k,4]],[g[k,3]],[g[k,5]]])
    savetxt('/Users/dmelgar/Slip_inv/Chiapas/data/statics/'+sta[k]+'.neu',neu,fmt='%.6f')