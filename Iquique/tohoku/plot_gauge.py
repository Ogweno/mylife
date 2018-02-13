from numpy import genfromtxt,arange,mean
from matplotlib import pyplot as plt

T=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/TI.MYE.txt')

t=arange(0,len(T)*15,15)/90.
tsun=T[:,6]
tsun=tsun-mean(tsun[1:10])

plt.figure()
plt.subplot(311)
plt.plot([25.5,25.5],[-10,10],'--',c='r')
plt.plot(t,tsun,'k')
plt.xlim([0,90])
plt.ylim([-0.75,1.0])
plt.legend(['Tsunami arrival'],frameon=False)
plt.ylabel('TI.MYE (m)')
plt.grid()



T=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/TI.OKA.txt')

t=arange(0,len(T)*15,15)/90.
tsun=T[:,6]
tsun=tsun-mean(tsun[1:10])

plt.subplot(312)
plt.plot(t,tsun,'k')
plt.plot([23.9,23.9],[-10,10],'--',c='r')
plt.xlim([0,90])
plt.ylim([-1,1.2])
plt.xlabel('Minutes after OT')
plt.ylabel('TI.OKA (m)')
plt.grid()

T=genfromtxt('/Users/dmelgar/Tohoku2011/Tsunami/TI.ONA.txt')

t=arange(0,len(T)*15,15)/90.
tsun=T[:,6]
tsun=tsun-tsun[0]

plt.subplot(313)
plt.plot(t,tsun,'b')
plt.plot([15.05,15.05],[-10,10],'--',c='r')
plt.xlim([0,90])
plt.ylim([-3.5,3.9])
plt.xlabel('Minutes after OT')
plt.ylabel('TI.ONA (m)')
plt.grid()
plt.show()