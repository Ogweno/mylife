from obspy import read
from matplotlib import pyplot as plt

stas=['coqu','pich','quin','valp']
path=u'/Users/dmelgar/Coquimbo2015/tsunami/sac/'
xl=[0,80]
yl=[-2.5,4.5]
xt=[15,30,45,60,75]

s0=read(path+stas[0]+'.sac')
s1=read(path+stas[1]+'.sac')
s2=read(path+stas[2]+'.sac')
s3=read(path+stas[3]+'.sac')

plt.figure(figsize=(3,7))

plt.subplot(411)
plt.plot(s0[0].times()/60,s0[0].data,'b',lw=1.5)
plt.grid()
plt.xlim(xl)
plt.ylim(yl)
plt.xticks(xt)
plt.yticks([-2,0,2,4])
plt.annotate('COQU',xy=(5,2.2))


plt.subplot(412)
plt.plot(s1[0].times()/60,s1[0].data,'b',lw=1.5)
plt.grid()
plt.xlim(xl)
plt.ylim(yl)
plt.xticks(xt)
plt.yticks([-2,0,2,4])
plt.annotate('PICH',xy=(5,2.2))

plt.subplot(413)
plt.plot(s2[0].times()/60,s2[0].data,'b',lw=1.5)
plt.grid()
plt.xlim(xl)
plt.ylim(yl)
plt.xticks(xt)
plt.yticks([-2,0,2,4])
plt.annotate('QUIN',xy=(5,2.2))

plt.subplot(414)
plt.plot(s3[0].times()/60,s3[0].data,'b',lw=1.5)
plt.grid()
plt.xlim(xl)
plt.ylim(yl)
plt.xticks(xt)
plt.yticks([-2,0,2,4])
plt.annotate('VALP',xy=(5,2.2))

plt.subplots_adjust(left=0.25,right=0.75,hspace=0.55)
plt.show()