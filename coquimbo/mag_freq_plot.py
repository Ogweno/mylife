from numpy import genfromtxt,diff,linspace,log10
import matplotlib.pyplot as plt


s=genfromtxt('/Users/dmelgar/Coquimbo2015/afters/csn_seismicity.txt')

n, bins, patches = plt.hist(s[:,6], 12)

delta=diff(bins)[0]

plt.figure()
plt.semilogy(bins[0:-1]+delta,n,'k')
plt.scatter(bins[0:-1]+delta,n,s=50,c='r')
#Plot reference line
a=log10(n[7])+(bins+delta)[7]
N=10**(a-linspace(1,8,100))
plt.semilogy(linspace(1,8,100),N,'--')
plt.grid(which='both')
plt.xlabel('Magnitude')
plt.ylabel('N')

plt.show()


